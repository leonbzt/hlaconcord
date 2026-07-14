"""Command-line interface (PLAN.md §10).

Thin orchestration over the importable pipeline: resolve a reference release,
build the in-house nomenclature backend, run parse -> normalize -> concord, and
write the tidy table / concordance report / JSON / GL strings.

    hlacc run --inputs optitype:S1.tsv arcashla:S1.json --db 3.55.0 -o out/
    hlacc run --samplesheet samples.csv -o out/
    hlacc validate  HLA-A*02:01 A*99:99
    hlacc normalize A*0201
    hlacc db list | update <ver> | pin <ver> | path
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import TextIO

from . import __version__, db, gl, report, tidy
from .concordance import ConsensusRule
from .nomenclature import (
    InHouseNomenclature,
    ParseError,
    ReductionBasis,
    ValidationStatus,
)
from .pipeline import InputSpec, PipelineResult, run

_BASIS_CHOICES = {b.value: b for b in ReductionBasis}
_RULE_CHOICES = {r.value: r for r in ConsensusRule}


# -- shared helpers -----------------------------------------------------------

def _root(args: argparse.Namespace) -> Path:
    return Path(args.data_dir).expanduser() if args.data_dir else db.default_root()


def _build_nomenclature(root: Path, requested: str | None) -> tuple[InHouseNomenclature, str]:
    """Resolve the release, load its reference data, return (backend, version)."""
    version = db.resolve_version(root, requested)
    reference = db.load_reference(root, version)
    return InHouseNomenclature(reference), version


def _fail(message: str, code: int = 2) -> int:
    print(f"hlaconcord: {message}", file=sys.stderr)
    return code


# -- `run` --------------------------------------------------------------------

def _load_samplesheet(path: Path) -> list[InputSpec]:
    """Read a batch samplesheet: columns sample,tool,path[,db_version]."""
    specs: list[InputSpec] = []
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        required = {"sample", "tool", "path"}
        missing = required - {(c or "").strip().lower() for c in (reader.fieldnames or [])}
        if missing:
            raise ValueError(
                f"samplesheet {path} missing column(s): {', '.join(sorted(missing))} "
                "(need sample,tool,path[,db_version])"
            )
        for row in reader:
            norm = {(k or "").strip().lower(): (v or "").strip() for k, v in row.items()}
            if not norm.get("tool") and not norm.get("path"):
                continue  # blank line
            base = path.parent
            file_path = Path(norm["path"])
            if not file_path.is_absolute():
                file_path = base / file_path
            specs.append(
                InputSpec(
                    tool=norm["tool"],
                    path=file_path,
                    sample=norm.get("sample") or None,
                    source_db_version=norm.get("db_version") or None,
                )
            )
    if not specs:
        raise ValueError(f"samplesheet {path} contains no input rows")
    return specs


def _collect_inputs(args: argparse.Namespace) -> list[InputSpec]:
    if args.samplesheet:
        return _load_samplesheet(Path(args.samplesheet))
    specs = [InputSpec.parse_arg(item) for item in args.inputs]
    if args.sample:
        for spec in specs:
            spec.sample = args.sample
    return specs


def _write_outputs(result: PipelineResult, out_dir: Path, *, want_gl: bool) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    gl_strings = gl.gl_strings_by_sample(result.concordance)
    meta = {
        "hlaconcord_version": __version__,
        "db_version": result.db_version,
        "basis": result.basis.value,
        "consensus_rule": result.rule.value,
        "tools": result.tools,
        "samples": result.samples,
    }
    written: list[Path] = []

    tidy_path = out_dir / "tidy.tsv"
    with tidy_path.open("w", newline="", encoding="utf-8") as fh:
        tidy.write_tidy_tsv(result.calls, fh)
    written.append(tidy_path)

    conc_path = out_dir / "concordance.tsv"
    with conc_path.open("w", newline="", encoding="utf-8") as fh:
        report.write_concordance_tsv(result.concordance, fh)
    written.append(conc_path)

    json_path = out_dir / "concordance.json"
    with json_path.open("w", encoding="utf-8") as fh:
        report.write_concordance_json(result.concordance, fh, meta=meta, gl_strings=gl_strings)
    written.append(json_path)

    if want_gl:
        gl_path = out_dir / "gl_strings.tsv"
        with gl_path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["sample_id", "gl_string"])
            for sample_id in result.samples:
                writer.writerow([sample_id, gl_strings.get(sample_id, "")])
        written.append(gl_path)

    return written


def _print_summary(result: PipelineResult, *, want_gl: bool, stream: TextIO) -> None:
    print(
        f"# hlaconcord {__version__}  db={result.db_version}  basis={result.basis.value}  "
        f"consensus={result.rule.value}  tools={','.join(result.tools)}",
        file=stream,
    )
    print(report.format_concordance(result.concordance), file=stream)
    if want_gl:
        print(file=stream)
        for sample_id, gl_string in gl.gl_strings_by_sample(result.concordance).items():
            print(f"# GL {sample_id}\n{gl_string or '(no consensus)'}", file=stream)


def _cmd_run(args: argparse.Namespace) -> int:
    root = _root(args)
    try:
        specs = _collect_inputs(args)
    except (ValueError, OSError) as exc:
        return _fail(str(exc))
    try:
        nom, version = _build_nomenclature(root, args.db)
    except db.DatabaseError as exc:
        return _fail(str(exc))

    try:
        result = run(
            specs,
            nom,
            db_version=version,
            basis=_BASIS_CHOICES[args.compare],
            rule=_RULE_CHOICES[args.consensus],
        )
    except (ValueError, OSError) as exc:
        return _fail(str(exc))

    if args.out:
        written = _write_outputs(result, Path(args.out), want_gl=args.gl)
        if not args.quiet:
            _print_summary(result, want_gl=args.gl, stream=sys.stdout)
            print("\n" + "\n".join(f"wrote {p}" for p in written))
    else:
        _print_summary(result, want_gl=args.gl, stream=sys.stdout)

    discordant = sum(1 for r in result.concordance if r.status.value == "discordant")
    return 1 if discordant else 0


# -- `validate` / `normalize` -------------------------------------------------

_OK_STATUSES = {
    ValidationStatus.EXACT,
    ValidationStatus.VALID_REDUCTION,
    ValidationStatus.G_GROUP,
    ValidationStatus.P_GROUP,
}


def _cmd_validate(args: argparse.Namespace) -> int:
    root = _root(args)
    try:
        nom, version = _build_nomenclature(root, args.db)
    except db.DatabaseError as exc:
        return _fail(str(exc))

    print(f"# db={version}")
    print("name\tstatus\taccession")
    exit_code = 0
    for name in args.names:
        try:
            allele = nom.parse(name)
        except ParseError as exc:
            print(f"{name}\tunparseable\t")
            print(f"  ! {exc}", file=sys.stderr)
            exit_code = 1
            continue
        status = nom.validate(allele)
        accession = nom.accession(allele) or ""
        print(f"{name}\t{status.value}\t{accession}")
        if status not in _OK_STATUSES:
            exit_code = 1
    return exit_code


def _cmd_normalize(args: argparse.Namespace) -> int:
    root = _root(args)
    try:
        nom, version = _build_nomenclature(root, args.db)
    except db.DatabaseError as exc:
        return _fail(str(exc))

    basis = _BASIS_CHOICES[args.compare]
    print(f"# db={version}")
    print("raw\tnormalized\tresolution\tnull\t" + basis.value + "\taccession")
    exit_code = 0
    for name in args.names:
        try:
            allele = nom.parse(name)
        except ParseError as exc:
            print(f"{name}\t(unparseable)\t\t\t\t")
            print(f"  ! {exc}", file=sys.stderr)
            exit_code = 1
            continue
        reduced = nom.reduce(allele, basis)
        accession = nom.accession(allele) or ""
        print(
            f"{name}\t{allele.name()}\t{allele.resolution.value}\t"
            f"{'yes' if allele.is_null else 'no'}\t{reduced}\t{accession}"
        )
    return exit_code


# -- `db` ---------------------------------------------------------------------

def _cmd_db_list(args: argparse.Namespace) -> int:
    root = _root(args)
    releases = db.list_releases(root)
    pinned = db.pinned_version(root)
    print(f"# data dir: {root}")
    if not releases:
        print("(no releases installed — run `hlacc db update <version>`)")
        return 0
    for version in releases:
        mark = "*" if pinned and db.version_of(pinned) == version else " "
        print(f"{mark} {version}")
    if not pinned:
        print("(none pinned)")
    return 0


def _cmd_db_path(args: argparse.Namespace) -> int:
    print(_root(args))
    return 0


def _cmd_db_update(args: argparse.Namespace) -> int:
    root = _root(args)
    try:
        result = db.update(
            root, args.version, with_history=not args.no_history, set_pin=not args.no_pin
        )
    except db.DatabaseError as exc:
        return _fail(str(exc))
    print(f"installed {result.version} -> {result.directory}")
    print(f"  files: {', '.join(result.files)}")
    if result.skipped_history:
        print("  note: skipped Allelelist_history.txt (no cross-version reconciliation)")
    if not args.no_pin:
        print(f"  pinned {result.version} as default")
    return 0


def _cmd_db_pin(args: argparse.Namespace) -> int:
    root = _root(args)
    try:
        pinned = db.pin(root, args.version)
    except db.DatabaseError as exc:
        return _fail(str(exc))
    print(f"pinned {pinned}")
    return 0


# -- parser -------------------------------------------------------------------

def _add_db_option(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--db", help="IPD-IMGT/HLA release to validate against (e.g. 3.55.0)")
    parser.add_argument("--data-dir", help="override the reference-database cache root")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="hlacc",
        description="Harmonize and validate HLA typing output across multiple typers.",
    )
    parser.add_argument("--version", action="version", version=f"hlaconcord {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    # run
    p_run = sub.add_parser("run", help="harmonize typer outputs into a concordance report")
    src = p_run.add_mutually_exclusive_group(required=True)
    src.add_argument("--inputs", nargs="+", metavar="TOOL:PATH", help="typer outputs to ingest")
    src.add_argument(
        "--samplesheet", metavar="CSV", help="batch: columns sample,tool,path[,db_version]"
    )
    p_run.add_argument("--sample", help="sample id (overrides names parsed from --inputs files)")
    p_run.add_argument(
        "--compare", choices=sorted(_BASIS_CHOICES), default=ReductionBasis.LGX.value,
        help="comparison basis (default: lgx)",
    )
    p_run.add_argument(
        "--consensus", choices=sorted(_RULE_CHOICES), default=ConsensusRule.MAJORITY.value,
        help="consensus rule (default: majority)",
    )
    p_run.add_argument("-o", "--out", help="output directory for tidy/concordance/json[/gl]")
    p_run.add_argument("--gl", action="store_true", help="also emit consensus GL strings")
    p_run.add_argument(
        "--quiet", action="store_true", help="suppress the stdout summary when writing files"
    )
    _add_db_option(p_run)
    p_run.set_defaults(func=_cmd_run)

    # validate
    p_val = sub.add_parser("validate", help="classify allele names against a release")
    p_val.add_argument("names", nargs="+", metavar="NAME")
    _add_db_option(p_val)
    p_val.set_defaults(func=_cmd_validate)

    # normalize
    p_norm = sub.add_parser("normalize", help="show canonical form + accession + reduction")
    p_norm.add_argument("names", nargs="+", metavar="NAME")
    p_norm.add_argument(
        "--compare", choices=sorted(_BASIS_CHOICES), default=ReductionBasis.LGX.value,
        help="reduction basis to display (default: lgx)",
    )
    _add_db_option(p_norm)
    p_norm.set_defaults(func=_cmd_normalize)

    # db
    p_db = sub.add_parser("db", help="manage local IPD-IMGT/HLA releases")
    db_sub = p_db.add_subparsers(dest="db_command", required=True)

    d_list = db_sub.add_parser("list", help="list installed releases")
    d_list.add_argument("--data-dir")
    d_list.set_defaults(func=_cmd_db_list)

    d_path = db_sub.add_parser("path", help="print the cache root")
    d_path.add_argument("--data-dir")
    d_path.set_defaults(func=_cmd_db_path)

    d_update = db_sub.add_parser("update", help="fetch a release from ANHIG/IMGTHLA")
    d_update.add_argument("version", help="release to fetch, e.g. 3.55.0")
    d_update.add_argument("--no-history", action="store_true", help="skip the large history file")
    d_update.add_argument("--no-pin", action="store_true", help="do not pin the fetched release")
    d_update.add_argument("--data-dir")
    d_update.set_defaults(func=_cmd_db_update)

    d_pin = db_sub.add_parser("pin", help="set the default release")
    d_pin.add_argument("version")
    d_pin.add_argument("--data-dir")
    d_pin.set_defaults(func=_cmd_db_pin)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    # `db` subcommands don't all take --db; default it so helpers can read args.db.
    if not hasattr(args, "db"):
        args.db = None
    if not hasattr(args, "data_dir"):
        args.data_dir = None
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
