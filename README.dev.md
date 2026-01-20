# `mibitrans` developer documentation

## Development install

```shell
# Create a virtual environment, e.g. with
python -m venv env

# activate virtual environment
source env/bin/activate

# make sure to have a recent version of pip and setuptools
python -m pip install --upgrade pip setuptools

# (from the project root directory)
# install mibitrans as an editable package
python -m pip install --no-cache-dir --editable .
# install development dependencies
python -m pip install --no-cache-dir --editable .[dev]
```

Afterwards check that the install directory is present in the `PATH` environment variable.

## Running the tests

There are two ways to run tests.

The first way requires an activated virtual environment with the development tools installed:

```shell
pytest -v
```

The second is to use `tox`, which can be installed separately (e.g. with `pip install tox`), i.e. not necessarily inside the virtual environment you use for installing `mibitrans`, but then builds the necessary virtual environments itself by simply running:

```shell
tox
```

Testing with `tox` allows for keeping the testing environment separate from your development environment.
The development environment will typically accumulate (old) packages during development that interfere with testing; this problem is avoided by testing with `tox`.

### Test coverage

In addition to just running the tests to see if they pass, they can be used for coverage statistics, i.e. to determine how much of the package's code is actually executed during tests.
In an activated virtual environment with the development tools installed, inside the package directory, run:

```shell
coverage run
```

This runs tests and stores the result in a `.coverage` file.
To see the results on the command line, run

```shell
coverage report
```

`coverage` can also generate output in HTML and other formats; see `coverage help` for more information.

## Running linters locally

For linting and sorting imports we will use [ruff](https://beta.ruff.rs/docs/). Running the linters requires an 
activated virtual environment with the development tools installed.

```shell
# linter
ruff check .

# linter with automatic fixing
ruff check . --fix
```

To fix readability of your code style you can use [yapf](https://github.com/google/yapf).

You can enable automatic linting with `ruff` on commit by enabling the git hook from `.githooks/pre-commit`, like so:

```shell
git config --local core.hooksPath .githooks
```

## Cleaning notebooks

The linter tests now also check that the example notebooks are "clean". This means they do not contain the output and execution metadata that
Jupyter adds when you execute them. Removing this reduces the size and noisiness of the diffs and consequently makes reviewing changes easier.

To clean the notebooks locally before each commit, you can use the `nb-clean` tool, which is listed as one of the `[dev]` dependencies of the project. They are installable
with:
```
python -m pip install .[dev]
```

You can then run:
```
nb-clean add-filter
```
This adds a git hook that will automatically clean any staged notebooks before they are committed. If you would rather run this manually, you can instead use:
```
nb-clean clean mynotebook.ipynb
```
replacing with the name of the notebook in question.




## Testing docs locally

To build the documentation locally, first make sure `mkdocs` and its dependencies are installed:
```shell
python -m pip install .[doc]
```

Then you can build the documentation and serve it locally with
```shell
mkdocs serve
```

This will return a URL (e.g. `http://127.0.0.1:8000/mibitrans/`) where the docs site can be viewed.

Note that this will only create the "non-versioned" documentation, which should be fine for testing changes to the docs.
The versioned documentation is created using the python utility called [mike](https://github.com/jimporter/mike?tab=readme-ov-file#mike) and its corresponding [mkdocs integration](https://squidfunk.github.io/mkdocs-material/setup/setting-up-versioning/).
In general it should not be necessary to test this, but if necessary, use [the mike documentation](https://github.com/jimporter/mike?tab=readme-ov-file#viewing-your-docs) to inspect locally.

## Versioning

Bumping the version across all files is done with [bump-my-version](https://github.com/callowayproject/bump-my-version), e.g.

```shell
bump-my-version bump major  # bumps from e.g. 0.3.2 to 1.0.0
bump-my-version bump minor  # bumps from e.g. 0.3.2 to 0.4.0
bump-my-version bump patch  # bumps from e.g. 0.3.2 to 0.3.3
```

## Making a release
To create a release you need write permission on the repository.

This section describes how to make a release:

1. preparation
1. making a release on GitHub

### (1/2) Preparation

1. Checkout the main branch locally
1. Verify that the information (especially the author list) in `CITATION.cff` is correct.
1. Make sure the [version has been updated](#versioning).
1. Run the unit tests with `pytest -v`
1. Make sure the [docs build and look good](#testing-docs-locally)


### (2/2) GitHub

When all is well, navigate to the [releases on GitHub](https://github.com/MiBiPreT/mibitrans/releases).

1. Press draft a new release button
1. Select the "Choose a tag" drop down and write out the new version (e.g. v1.3.2)
1. Press "Generate release notes" to automatically fill the title (with the version number) and generate a description (the changelog from the merge pull requests)
1. Press the Publish release button

This will create the release on github and automatically trigger:

1. The `.github/workflows/publish.yml` workflow which will build the package and publish it on PyPI
1. The Zenodo-Github integration into making a snapshot of your repository and sticking a DOI on it and adding the new version to the main Zenodo entry for your software.