[![github repo badge](https://img.shields.io/badge/github-repo-000.svg?logo=github&labelColor=gray&color=blue)](https://github.com/MiBiPreT/mibitrans)[![PyPI version](https://badge.fury.io/py/mibitrans.svg)](https://badge.fury.io/py/mibitrans)[![github license badge](https://img.shields.io/github/license/MiBiPreT/mibitrans)](https://github.com/MiBiPreT/mibitrans) [![RSD](https://img.shields.io/badge/rsd-mibitrans-00a3e3.svg)](https://www.research-software.nl/software/mibitrans) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10877886.svg)](https://doi.org/10.5281/zenodo.10877886) [![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow)](https://fair-software.eu) [![workflow scq badge](https://sonarcloud.io/api/project_badges/measure?project=MiBiPreT_anatrans&metric=alert_status)](https://sonarcloud.io/dashboard?id=MiBiPreT_anatrans) [![workflow scc badge](https://sonarcloud.io/api/project_badges/measure?project=MiBiPreT_anatrans&metric=coverage)](https://sonarcloud.io/dashboard?id=MiBiPreT_anatrans) [![documentation](https://github.com/MiBiPreT/mibitrans/actions/workflows/documentation-deploy-release.yml/badge.svg)](https://mibipret.github.io/mibitrans) [![build](https://github.com/MiBiPreT/mibitrans/actions/workflows/build.yml/badge.svg)](https://github.com/MiBiPreT/mibitrans/actions/workflows/build.yml) [![cffconvert](https://github.com/MiBiPreT/mibitrans/actions/workflows/cffconvert.yml/badge.svg)](https://github.com/MiBiPreT/mibitrans/actions/workflows/cffconvert.yml) [![sonarcloud](https://github.com/MiBiPreT/mibitrans/actions/workflows/sonarcloud.yml/badge.svg)](https://github.com/MiBiPreT/mibitrans/actions/workflows/sonarcloud.yml) [![markdown-link-check](https://github.com/MiBiPreT/mibitrans/actions/workflows/markdown-link-check.yml/badge.svg)](https://github.com/MiBiPreT/mibitrans/actions/workflows/markdown-link-check.yml)

`mibitrans` is an analytical subsurface contaminant transport modelling tool based on the three-dimensional advection–dispersion equation. `mibitrans` has a clear and quick parameter input and various model options, which allow for easy visualization of contaminant transport under various field conditions. The package is designed to be modular and well documented, to allow users to adapt a model to their liking. `mibitrans` is develloped as a screening tool to gain insights into field site contaminant distribution, which then can be used in development of more involved numerical models. The package is also well suited in contaminant transport visualization for educational purposes.

## Installation of stable release

Use `pip` to install the most recent stable release of `mibitrans` as follows:

```console
pip install mibitrans
```

## Installation of most recent development version

To install mibitrans from the GitHub repository directly, do:

```console
git clone git@github.com:MiBiPreT/mibitrans.git
cd mibitrans
python -m pip install .
```

Note that this is the (possibly unstable) development version from the `main` branch. If you want a stable release, use the pip installation method instead.

## Documentation

[See the full `mibitrans` documentation here](https://mibipret.github.io/mibitrans/)

## Contributing

If you want to contribute to the development of mibitrans,
have a look at the [contribution guidelines](CONTRIBUTING.md).

## Credits

This package was created with [Cookiecutter](https://github.com/audreyr/cookiecutter) and the [NLeSC/python-template](https://github.com/NLeSC/python-template).
