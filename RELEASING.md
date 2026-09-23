# Releasing pyoma

pyoma's version is *not* stored in a file — it is derived from git tags via
[`setuptools_scm`](https://setuptools-scm.readthedocs.io/) (see
`[tool.setuptools_scm]` in `pyproject.toml`). There is nothing to bump by hand:

- A checkout exactly on tag `v1.2.3` builds as `1.2.3`.
- Any other commit builds as something like `1.2.4.dev6+gabc1234`
  (next patch, 6 commits past the last tag, at commit `abc1234`), which is
  PEP 440-compliant and always sortable/installable, but not something
  you'd knowingly publish to PyPI.

This means the *only* releasable commits are tagged commits on `master`.

## Branch roles

- `develop` — integration branch. Feature branches merge here.
- `master` — release branch. Only receives merges from `develop` and release
  tags. Every push of a `v*` tag triggers
  `.github/workflows/upload-release.yml`, which builds and publishes to PyPI.

## Cutting a release

1. Make sure `develop` is green and has everything you want released.
2. Merge `develop` into `master`:
   ```
   git checkout master
   git pull
   git merge develop
   ```
3. Tag the merge commit with the new version (`v` prefix, semver):
   ```
   git tag vX.Y.Z
   git push origin master vX.Y.Z
   ```
   Pushing the tag triggers the PyPI upload workflow — no manual version
   edit, no separate "bump version" commit.
4. **Merge `master` back into `develop` immediately**, so the tag is part of
   `develop`'s history too:
   ```
   git checkout develop
   git merge master
   git push origin develop
   ```
   This step is not optional — skipping it is what caused `develop` and
   `master` to drift apart in the past (`develop` ended up missing several
   `master`-only commits, and `setuptools_scm` on `develop` would compute dev
   versions against a stale tag). Do this step every time, right after
   tagging.

## Choosing the version number

Standard semver against the previous tag:
- **patch** (`0.14.2` → `0.14.3`): bug fixes only.
- **minor** (`0.14.2` → `0.15.0`): new features, backwards compatible.
- **major** (`0.14.2` → `1.0.0`): breaking changes.

## Checking what a build will produce

Before tagging, you can preview the version any commit would build as:
```
uv run --with setuptools_scm python -m setuptools_scm
```
