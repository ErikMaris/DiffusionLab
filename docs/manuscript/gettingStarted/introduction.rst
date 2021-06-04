Introduction
=================

In tracking experiments -- common in the field of chemistry, biology, and physics -- time-lapse movies are obtained of to-be-tracked object(s) of interest. This is often molecules, proteins, or colloids. Dedicated software can record the location of the object in each frame and group the positions over time that belong to the same object. The groups of locations are called tracks (or trajectories) and contain information about, among others, the motion of these molecules or particles. This motion can contain value information about the tracked object itself and its surroundings, but this is often not straightforward to obtain for either two reasons:

1. Trajectories contain only a few locations due to fast motion out of focus or limited stability of the tracked object, and can therefore not be analyzed individually.
2. The trajectories with different fundamental behavior are present in the same data set.

**DiffusionLab** is software that was designed to deal with the problems described above, but has grown into a versatile track motion analysis software. By first performing a segmentation of the tracks into user-defined categories, erroneous tracks can be removed and tracks with the same motion behavior can be grouped. The model can be tailored to each group individually and computed for the tracks. The mean motion behavior of a group can be computed if the tracks are too short to be analyzed individually. The suggested workflow can be found :doc:`here <workflow>`.