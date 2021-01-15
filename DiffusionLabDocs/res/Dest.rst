.. _ch-diffusionEstimators:

Diffusion estimators
------------------------

In this resource section, the available diffusion estimators are documented.

The motion blur coefficient :math:`R` is computed following equation S9 in Lindén et al. :cite:`linden2018variational`. A continuous exposure at the beginning of each frame is assumed.

MSD
++++

Mean-squared displacement (MSD) analysis is performed both on the individual tracks MSD as the population MSD. As the number of displacements contributing to the mean decreases at longer delay times and the bias increases, it is recommended to not fit the full MSD curve. The fit range can be set in the :guilabel:`Diffusion estimator options` menu.

* **Clipping factor**: maximum fraction of the total number of delay times, when value < 1 or the number of delay times included in the fit when value > 1. Please note that a delay time of zero is excluded from the fit.
* **Minimum points taken for fit**: minimum number of delay times taken for fit regardless of the clipping factor.

Diffusion-estimator specific plots:

* **Plot MSD fit of track**: plots the MSD fit of the current track.
* **plot MSD fit of all tracks**: plots the population MSD fit of the current population.


MSD_normal
==============

Model for normal Brownian motion. Fits the fit range of the MSD :math:`<r^2>` with :cite:`michalet2012optimal`:

.. math::
	<r^2> = 2dDt + 2d\sigma^2 - 4dRD\Delta t

* :math:`d`: dimension
* :math:`D`: diffusion constant
* :math:`t`: delay time
* :math:`\sigma`: localization error
* :math:`R`: motion blur constant
* :math:`\Delta t`: frame time

The diffusion signal-to-noise (SNR) is calculated following Vestergaard et al. :cite:`vestergaard2014optimal`:

.. math:
	\textrm{SNR} = \frac{\sqrt{D \Delta t}}{\sigma}

The :math:`D`, :math:`\sigma^2`, and :math:`\textrm{SNR}` are added as track properties.

MSD_confined
==============

Model for trapping in domains. Fits the fit range of the MSD :math:`<r^2>` with :cite:`saxton2007modeling`:

.. math::
	<r^2> = <r^2>_0 \left[ 1 - \exp(-t/\tau) \right[

* :math:`<r^2>_0`: squared confinement length
* :math:`t`: delay time
* :math:`\tau`: confinement time

The :math:`<r^2>_0` and :math:`\tau` are added as track properties. More information on the relation between the squared confiment length and the confiment domain size can be found in Qian el al. :cite:`qian1991single`.

MSD_directed
==============

Model for normal Brownian motion in a flow. Fits the fit range of the MSD :math:`<r^2>` with :cite:`saxton2007modeling`:

.. math::
	<r^2> = 2dDt + 2d\sigma^2 - 4dRD\Delta t + v^2t^2

* :math:`d`: dimension
* :math:`D`: diffusion constant
* :math:`t`: delay time
* :math:`\sigma`: localization error
* :math:`R`: motion blur constant
* :math:`\Delta t`: frame time
* :math:`v`: velocity

The :math:`D`, :math:`\sigma^2`, :math:`v`, and :math:`\textrm{SNR}` are added as track properties.

