# DiffusionLab
Software package for single-fluorophore trajectory analysis.

This code was inspired by [MSD analyzer]{https://tinevez.github.io/msdanalyzer/}.

Link to documentation: https://diffusionlab.readthedocs.io/en/latest/.

## Requirements

Software has been written and tested in MATLAB 2019b.

(At least) the following MATLAB toolboxes must be installed:
- Curve Fitting Toolbox
- Image Fitting Toolbox

When running from the command line, please make sure that the following repositories have been downloaded and added to the MATLAB path:
- printFig (https://github.com/ErikMaris/printFig)
- plotProps (https://github.com/ErikMaris/plotProps)
- unitProps (https://github.com/ErikMaris/unitProps)

## Demo data

The trajectories from *Maris, J.J.E. et al (in preparation) 2021* are available to test the software installation. The trajectories are available in  the folder 
```
data/processed/DiffusionLab_manuscript/
```
with the name 
```
example_1-2_tracks_sig2_SNR3_FP1_BG2_PJ20.csv
```
and a training set for the Classification Trainer (for DiffusionLab) can be found in the same folder with the name 
```
1-training.mat
```
. The trajectories in 
```
example_1-1_trainingset_tracks_sig2_SNR3_FP1_BG2_PJ15.csv
```
were used to construct the training set.

Please open the DiffusionLab App and load 
```
example_1-2_tracks_sig2_SNR3_FP1_BG2_PJ20.csv
```
in DiffusionLab via *File* > *Import tracks* > *DoM*: this should result in 10600 imported trajectories. First set the *Pixel size (nm)* to 64, *Frame time (s)* 0.050, and *Exposure time (s)* to 0. The shortest trajectories should be removed first to facilitate further analysis. Select *Number of Points* as first track property (below *Track 1*) and set a *Property threshold* of 5. Press *Use property 1 as filter* to select the threshold and *Apply filter to all tracks* to execute. Delete the first population containing the tracks with fewer than 5 points by pressing *Delete current population*. The remain population should have 2433 trajectories. Now, the data set can be visualized, for instance by plotting the tracks using *Plot tracks* or their mean squared displacement curves with *Plot MSD all tracks*.

For more information about the usage of the DiffusionLab graphical user interface please consult the [documentation](https://diffusionlab.readthedocs.io/en/latest/).

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
