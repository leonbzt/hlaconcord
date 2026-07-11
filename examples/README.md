# Example: harmonizing one sample across four typers

`s1/` holds the raw output of the four MVP typers for a single sample, each in
its own format and resolution — the exact heterogeneity hlaharm exists to absorb:

| File | Tool | Resolution emitted |
|---|---|---|
| `s1_result.tsv` | OptiType | 2-field (`A*02:01`) |
| `s1.genotype.json` | arcasHLA | 3-field (`A*02:01:01`) |
| `s1_bestguess_G.txt` | HLA-LA | G-group (`A*02:01:01G`) |
| `s1_final.result.txt` | HLA-HD | 3-field, `HLA-` prefixed (`HLA-A*02:01:01`) |

## Run it

First install a reference release (fetched into a local cache, not bundled):

```bash
hlaharm db update 3.55.0
```

Then harmonize the four outputs:

```bash
hlaharm run \
  --inputs optitype:s1/s1_result.tsv \
           arcasHLA:s1/s1.genotype.json \
           hla-la:s1/s1_bestguess_G.txt \
           hla-hd:s1/s1_final.result.txt \
  --db 3.55.0 --gl -o out/
```

## Expected result

All four resolutions of each allele collapse to one `lgx` key, and every locus is
concordant:

```
# sample s1 (basis=lgx)
  OK A      A*01:01 + A*02:01
  OK B      B*07:02
  OK C      C*07:02
  OK DRB1   DRB1*15:01

# GL s1
HLA-A*01:01+HLA-A*02:01^HLA-B*07:02^HLA-C*07:02^HLA-DRB1*15:01
```

`out/` then contains `tidy.tsv` (one auditable row per allele call, raw alongside
normalized), `concordance.tsv`, `concordance.json`, and `gl_strings.tsv`.

DRB1 is called by three tools — OptiType is class I only, so its absence there is
"not called", not a disagreement.
