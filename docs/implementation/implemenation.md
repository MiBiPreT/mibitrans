# Implemenation

## Architecture

The code is structured in a modular fashion to ensure easy extension of functionality by the (future) development team as well as outside collaborators. To stimulate embedding in the scientific community, users are encouraged to report bugs and suggest extensions and issues on the GitHub issue tracker, or even self-implemented features, through pull requests.

<img src="mibitrans_schematics_structure.png" alt="Overview of `mibitrans` package structure and functionalities." width="500">

The general structure follows the split between input (data module), process control (transport module), and output (visualize module). Input data follows in sorting the related processes: hydrological/aquifer parameters, attenuation parameters, source parameters, and model specific parameters. 

Transport of contaminants is solved analytically with a choice of 3 model options, corresponding to the solutions outlined in section **Background/Model implementations**. Each model is implemented in a separate class:

* `Mibitrans`: this model class provides the concentration distribution according to the fully exact solution. This solution requires numerical integration which causes slightly higher computational effort and duration of calculation but provides the most accurate results. Note the deliberate choice of naming the model class for the exact solution `Mibitrans` similar to the python package name `mibitrans` as we consider this model to be the default option.     
* `Anatrans: this model class provides the concentration distribution according to the analytical approximation. This makes numerical integration no longer necessary. However, it introduces an error whose size dependings on field conditions. The purpose of this model is to fasten up repeated model use, e.g. when used for parameter calibration.
* `Bioscreen: this model class provides the concentration distribution according to the approximate equation implemented in *BIOSCREEN*. The purpose of this model class is mostly benchmarking and comparison.

Transport models are implemented as classes, each inheriting base functionality, such as data input and results handling. This setup allows to simply define other transport model solutions, as we show in the discussion section. 

Currently, two options for modeling contaminant degradation are available: (i) linear decay and (ii) instant reaction (section **Background/Instant Reactin Model**). Standard degradation process is linear decay. The case of no decay is naturally included in this model implementation using a zero decay rate or half life time. The handling of degradation options is also implemented at transport parent class level. Both degradation options are thus available for each of the implemented transport classes.

A broad range of visualization options is implemented within the package. This does not only include display of calculated concentrations in various forms, but also visualization of the source settings. Examples will be provided for each visualization option in the following sections.

The `mibitrans` package provides the feature of calculating the mass balances within the source and mass within the plume. For instance, the remaining mass in the plume and the degraded mass can be calculated at each time.

In the section **examples**, various example workflows, including application examples, provide inside into the application of `mibitrans`-functionalityin in the form of commented Jupyter-Notebooks.

## Research Software Management

We wish to reach with `mibitrans` a wide range of academic researchers with their research and education as well as societal stakeholders, like remediation consultants, to quickly evaluate contaminated field sites. To achieve that we follow in the development of `mibitrans`, the three aspects:

* scientific usefulness
* FAIR principles for Research Software
* software quality

`mibitrans` is written in Python, to align with the predominant trend of using Python for geoscientific modeling and data analysis. To further instigate adoption by the scientific community and other stakeholders, the development team follows the latest standards for the development and management of research software.

The open-source software package is being developed on a public GitHub repository, licensed under the permissive Apache 2.0 license. It applies version control, with released versions published on Zenodo and the Python Package Index (PyPI). 
This last publication type enables a single-line installation with the \texttt{pip} tool, resulting in a very low usage threshold. 

To ensure high quality and consistent code style and scientific correctness, new additions are only integrated after passing automated unit tests and code quality checks and after passing a code review by a second research software engineer. 
The package has extensive documentation, for both users and developers, which is updated and published with every (minor) software release. 

