Conda environments used in the analyses.

[notebook_env.yml](./notebook_env.yml) is basically the [Best practices tutorial](https://github.com/theislab/single-cell-tutorial) environment with [scVelo](https://pypi.org/project/scvelo/) added. Check the tutorial [GitHub](https://github.com/theislab/single-cell-tutorial) for better installation instructions, including installation of some required R packages. This environment is used with the analysis notebooks.

[scanpy_env.yml](./scanpy_env.yml) includes almost the same packages as [notebook_env.yml](./notebook_env.yml) but with [mnnpy](https://pypi.org/project/mnnpy/) added. Also, originally scanpy was installed through GitHub [(version 1.4)](https://github.com/theislab/scanpy/tree/0401341877e586e36412db1e178601bf1f545037) with pip separately and is therefore not included in this environment file. See [here](https://github.com/conda/conda/issues/6805) for possible problems if wanting to include the automatic github install through pip in the conda environment. This environment is used when running the separate scanpy analyses on command line.

[pyscenic_env.yml](./pyscenic_env.yml) includes the packages needed to run the SCENIC analyses.

For maximum reproducibility, the conda environments are exported with full version tags. In different systems, the version tags have to be excluded from the environment files before installation.