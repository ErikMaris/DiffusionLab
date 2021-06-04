Export for DiffusionLab
========================

DiffusionLab supports the import of tracks from DoM v.1.1.6 and Localizer for Igor Pro. Support requests for formats can be sent to the developers.

DoM v.1.1.6
------------

DoM is an ImageJ plugin for localization and tracking of single-molecule fluorescence time lapse movies. The software is freely available `here <https://github.com/ekatrukha/DoM_Utrecht/wiki>`_. The results can be exported via the ImageJ “Results” table. The steps are the following:

1. Open your movie in ImageJ
2. Analyze > DoM v.1.1.6 > Detect molecules and perform “Detect molecules and fit” with the optimized settings for your system. Use a “fitting iterations number” of 5
3. Analyze > DoM v.1.1.6 > Link Particles to Tracks
4. Go to the “Results” table and go to File > Save As and save the table as .csv file

Macros are available in the group for batch localization and linking into tracks.


Localizer for Igor Pro
------------------------

`Localizer <https://doi.org/10.1117/1.JBO.17.12.126008>`_ is a plugin for Igor Pro for localization and tracking of single-molecule fluorescence time lapse movies. The software is freely `available here <https://bitbucket.org/pdedecker/localizer/src/master/>`_. Tracks can be computed and exported following:

1. Localizer > Read CCD data > Read CCD data from disk...
2. Press >> on the movie viewer and perform Localization and Tracking
3. Localizer > Manipulate Particle Tracks > Save tracks to text file...

Supported Localization algorithms are:
 
* Gaussian fitting
* Gaussian fitting (fixed width)
* Ellipsoidal Gaussian fitting
* MLEwG

.. warning::
	The localization algorithm is recognized by the number of columns in the output file. If this changes, the columns are interpreted incorrectly.

COMSOL
--------

`COMSOL <https://www.comsol.nl/>`_ is a commercial software packing for multiphysics modelling and allows single-particle trajectory simulation. For this, the `Particle tracing module <https://www.comsol.nl/particle-tracing-module>`_ should be installed. Tracks can be exported follwing:

1. To do!


.. warning::
	COMSOL exports the coordinates in three dimensions. The functionality of DiffusionLab for 3D trajectories is limited and some track properties are only computed in the XY-plane. A warning in the MATLAB log window appears when this occurs, and explains which track properties have been affected.