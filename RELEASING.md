# Releasing hlaconcord

Releases are published to PyPI by `.github/workflows/release.yml`, triggered by a
version tag (`v*`). Publishing uses **PyPI Trusted Publishing** (OIDC) — there is
no API token stored anywhere.

## One-time setup (do this once, before the first release)

1. **Create a PyPI account** at <https://pypi.org> (if you don't have one) with 2FA.

2. **Register a pending trusted publisher** on PyPI. Because the project doesn't
   exist on PyPI yet, add it as a *pending* publisher:
   PyPI → *Your account* → *Publishing* → *Add a new pending publisher* → fill in:

   | Field | Value |
   |---|---|
   | PyPI Project Name | `hlaconcord` |
   | Owner | `leonbzt` |
   | Repository name | `hlaconcord` |
   | Workflow name | `release.yml` |
   | Environment name | `pypi` |

   The **environment name must be exactly `pypi`** — it matches `environment: pypi`
   in `release.yml`.

3. **Create the `pypi` environment** on GitHub:
   repo → *Settings* → *Environments* → *New environment* → name it `pypi`.
   Optionally add a required reviewer so a release must be approved before it publishes.

## Cutting a release

1. Bump the version — edit `__version__` in `src/hlaconcord/__init__.py` (single
   source of truth; `pyproject.toml` reads it dynamically). Follow SemVer.
2. Update `CHANGELOG.md` (move items under a new `## [x.y.z] — DATE` heading).
3. Commit on `main` and make sure CI is green.
4. Tag and push:

   ```bash
   git tag v0.1.0
   git push origin v0.1.0
   ```

   The tag triggers `release.yml`, which builds the sdist + wheel, checks the
   metadata, verifies the tag matches the package version, and publishes to PyPI.
   (If you added a required reviewer to the `pypi` environment, approve the run.)

5. Optionally create a GitHub Release from the tag:

   ```bash
   gh release create v0.1.0 --generate-notes
   ```

## Optional: dry-run on TestPyPI first

To rehearse without touching real PyPI, register a second pending publisher on
<https://test.pypi.org> (same fields), add `repository-url:
https://test.pypi.org/legacy/` to the `pypa/gh-action-pypi-publish` step in a
throwaway branch/workflow, and push a pre-release tag. Then `pip install
--index-url https://test.pypi.org/simple/ hlaconcord` in a clean venv and run the
`examples/` walkthrough.
