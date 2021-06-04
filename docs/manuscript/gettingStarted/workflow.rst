Workflow
==========

The suggested workflow is the following:

1. Localization, tracking, and export in DoM or Localizer
2. Import in DiffusionLab
3. Compute track properties
4. Plot track(s) and/or properties and identify desired populations
5. Segment in desired populations and use the plot previews to assess the segmentation
6. Compute the diffusion for all populations with the most desired diffusion estimator
7. Plot results and save figures

When the trajectories contain too few locations to be analyzed individually with a reasonable error, it is recommended to analyze the mean motion behavior of the population. Both the mean squared displacement (MSD) and cumulative probability density (CPD) analysis allow the fitting of the population mean.

Any steps 3--7 can be skipped if only partial functionality is needed.