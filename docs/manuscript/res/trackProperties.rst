.. _ch-trackProperties:

Track Properties
=================

Track properties are features that describe a property of a trajectory. This can relate to e.g. the shape of the trajectory or number of points. In DiffusionLab, these track properties are used to segment the tracks into smaller populations with similar motion behavior.

.. note::
	A good set of track properties for segmentation describes properties that are different between the intended populations. DiffusionLab provides a large selection of track properties, and depending on the intended populations the best subset of track properties can differ between data sets.


.. _ch-trackProperties-standard:

Standard Track Properties
--------------------------

The standard track properties are computed in the same way regardless the settings in DiffusionLab via :guilabel:`Compute properties`.

Number of points
++++++++++++++++++++++

**Description:** number of consecutive localizations.

**Physical interpretation:** 

* Mobility
* Photostability
* Statistical significance

**Units:** -.

Length
++++++++++++++++++++++

**Description:** the sum of all the individual displacements.

**Physical interpretation:** 

* Mobility
* Photostability

**Units:** length.

MinBoundCircleRadius
++++++++++++++++++++++

**Description:** radius of the smallest enclosing circle that can be drawn around the localization coordinates, i.e. minimum bounding circle (MBC).

**Physical interpretation:** 

* spatial extension of localization coordinates

**Units:** length.


MBC minus CoM
++++++++++++++++++++++

**Description:** the distance between the center of the ``MinBoundCircleRadius`` and the center of mass, calculated as a percentage of the MBC radius. [#f1]_

**Physical interpretation:** 

* Evenness spatial distribution localizations. It gives an indication of how homogeneously points are spatially distributed.

**Units:** -.


Entropy
++++++++++++++++++++++

**Description:** Shannon’s entropy of the distribution of the localization coordinates within the enclosing square that is defined by two times ``MinBoundCircleRadius``. [#f2]_

**Physical interpretation:** 

* Statistical measurement of spatial randomness

**Units:** -.


Tortuosity
++++++++++++++++++++++

**Description:** the ratio of the distance between start and end points versus the length of the track.

**Physical interpretation:** 

* start-to-end directionality.

**Units:** -.

Elongation
++++++++++++++++++++++

**Description:** weight of the first principal component of localization coordinates.

**Physical interpretation:** 

* Directionality of localizations

**Units:** -.

Elongation angle
++++++++++++++++++++++

**Description:** direction of the first principal component of localization coordinates

**Physical interpretation:** 

* Direction of localizations

**Units:** -.

Other Track Properties
---------------------------

The optional track properties are specified in the diffusion constant estimator and are dependent on the settings thereof. 

Diffusion constant
+++++++++++++++++++++++++++

**Description:** magnitude of the diffusion.

**Units:** length^2/time.


Localization error
++++++++++++++++++++++

**Description:** imprecision in the localization. The deviation of a localization  estimate from its true position is ideally normally distributed in one dimension. The localization error is defined as the standard deviation of this normal distribution.


**Units:** length.

Diffusion SNR
++++++++++++++++++++++

**Description:** signal-to-noise (SNR) of the displacements as given in Vestergaard et al. [#f3]_

**Physical interpretation:** 

* relative magnitude of diffusion to the localization error

**Units:** -.


Underlying Descriptors
----------------------------

The standard track properties categorized by their main descriptors are given in :numref:`Table  %s <tab-underlying-descr>`.

.. _tab-underlying-descr:

.. list-table:: Standard track properties categorized by their main underlying descriptor.
   :widths: 25 50
   :header-rows: 1

   * - Descriptor
     - Track property
   * - Mobility, photostability
     - Number of points, length
   * - Spatial directionality
     - Tortuosity, elongation, elongation angle
   * - Uniformity spatial distribution
     - Minimum bounding circle radius, MBCC minus CoM, entropy
	 
.. rubric:: References (in footnotes)

.. [#f1] Hendriks, F.C., Meirer, F., Kubarev, A.V., Ristanović, Z., Roeffaers, M.B., Vogt, E.T., Bruijnincx, P.C. and Weckhuysen, B.M., 2017. Single-molecule fluorescence microscopy reveals local diffusion coefficients in the pore network of an individual catalyst particle. Journal of the American Chemical Society, 139, pp.13632-13635.
.. [#f2] Same as ref. 1.
.. [#f3] Vestergaard, C.L., Blainey, P.C. and Flyvbjerg, H., 2014. Optimal estimation of diffusion coefficients from single-particle trajectories. Physical Review E, 89, p.022726.