# Publishing FastNSA

## Account Setup

The PyPI Trusted Publisher must identify owner `fastnsa`, repository
`FastNSA`, workflow filename `publish.yml`, and environment `pypi`.
In GitHub's `pypi` environment, configure a **Tag** rule named `v*`, not a
Branch rule. Zero matching tags is normal before the first release.

Add Required reviewers if upload should require human approval. A tag rule
alone is not an approval gate. A sole maintainer must not enable Prevent
self-review unless another eligible reviewer is available. No PyPI API token
secret is required.

## Dry Run

After the workflow is committed and pushed, open **Actions > Prepare and
publish FastNSA > Run workflow** on `main`.

Every manual dispatch is a dry run, even if a tag is selected. It builds and
tests the platform wheels, builds and smoke-tests the source archive, checks
distribution metadata, and assembles a `fastnsa-release-0.3.0` artifact.
It never enters the `pypi` environment, uploads a package, or creates a GitHub
Release. The upload jobs are expected to be skipped.

The artifact contains:

- `dist/`: exactly four platform wheels and one source distribution.
- `SHA256SUMS`: checksums of the distribution files.
- `manifest.json`: version, source commit/ref, target platforms, and hashes.
- `RELEASE_NOTES.md`: the version-specific release notes.

The local wheel used to smoke-test the source archive is not uploaded.
Validation rejects universal `any` wheels and unrepaired `linux_*` wheels.

## Release

1. Confirm that the default-branch CI and a dry run both passed.
2. Ensure package metadata, CMake, source fallback version, and release notes
   agree. The current target is `0.3.0`; create no tag for a failed dry run.
3. Create and push `v0.3.0` at the intended commit on `main`. A tag must point
   to a commit reachable from `main`; do not move an existing release tag.
4. The tag push starts a fresh build/test run at that commit. A mismatched tag,
   missing platform, or failed check prevents publishing.
5. Inspect the assembled artifact, then approve the `pypi` environment if
   Required reviewers is configured. Without that protection, successful tag
   runs proceed directly to uploading.
6. The upload job verifies the artifact handoff and checks that the version
   does not already exist on PyPI. Only this job receives OIDC permission.
7. After successful PyPI upload, the workflow creates a **draft** GitHub
   Release with the same artifacts and checksums. Inspect and publish the draft
   manually when ready to announce the release.
8. In a clean environment, install the exact version from PyPI and run
   `python -m fastnsa._selftest`. Confirm all expected files are on PyPI before
   enabling the QGIS candidate workflow.

Do not merely publish a GitHub Release and expect a package upload: this
workflow is triggered by a matching tag push, not the Release event.

## Failure and Recovery

The workflow does not use `skip-existing`, overwrite PyPI files, delete
releases, or move tags. Any existing PyPI version stops the upload job,
including a partially uploaded release. Stop and compare the existing PyPI
file hashes against the approved artifact before deciding how to recover.
Do not blindly rerun an upload with different files or recreate the tag.

If a draft Release was partly created, inspect it before retrying the GitHub
Release job; the workflow deliberately does not overwrite an existing draft.
Checksums verify integrity, not independent source authenticity.

## Local Checks

Release-tool tests require `packaging` and `PyYAML`, not the native library:

```bash
python -m pip install packaging PyYAML
python -m unittest discover -s tests/release -t . -v
python tools/ci/release.py version
```

The `version`, `assemble`, `verify`, and `check-index` subcommands never
upload or create tags/releases. `check-index` makes a read-only PyPI request.
Actual GitHub Actions, OIDC, environment approval, and four-platform validation
must be checked on GitHub; local tests cannot certify account settings.

References: [PyPI Trusted Publishing](https://docs.pypi.org/trusted-publishers/using-a-publisher/)
and [GitHub environment rules](https://docs.github.com/en/actions/how-tos/deploy/configure-and-manage-deployments/manage-environments).
