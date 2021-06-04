# DiffusionLab
Software package for single-fluorophore trajectory analysis

Link to documentation: https://diffusionlab.readthedocs.io/en/latest/

## Requirements

This code was written and tested in MATLAB 2019b.

(At least) the following MATLAB toolboxes must be installed:
- Curve Fitting Toolbox
- Image Fitting Toolbox

When running from the command line, please make sure that the following repositories have been donwloaded:
- printFig (https://github.com/ErikMaris/printFig)
- plotProps (https://github.com/ErikMaris/plotProps)
- unitProps (https://github.com/ErikMaris/unitProps)


## Project organization
- PG = project-generated
- HW = human-writable
- RO = read only
```
.
├── .gitignore
├── CITATION.md
├── LICENSE.md
├── README.md
├── requirements.txt
├── bin                <- Compiled and external code, ignored by git (PG)
│   └── external       <- Any external source code, ignored by git (RO)
├── config             <- Configuration files (HW)
├── data               <- All project data, ignored by git
│   ├── processed      <- The final, canonical data sets for modeling. (PG)
│   ├── raw            <- The original, immutable data dump. (RO)
│   └── temp           <- Intermediate data that has been transformed. (PG)
├── docs               <- Documentation notebook for users (HW)
│   ├── manuscript     <- Manuscript source, e.g., LaTeX, Markdown, etc. (HW)
│   └── reports        <- Other project reports and notebooks (e.g. Jupyter, .Rmd) (HW)
├── results
│   ├── figures        <- Figures for the manuscript or reports (PG)
│   └── output         <- Other output for the manuscript or reports (PG)
└── src                <- Source code for this project (HW)

```


## License

This project is licensed under the terms of the [GPL-3.0 License](/LICENSE.md)
