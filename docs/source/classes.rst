Classes
=========
The :mod:`src` module contains all QUPS `classes`. 

* :class:`ChannelData` : represents the echo signal data and its time axis.
* :class:`UltrasoundSystem` : represent the combination of a :class:`Transducer`, :class:`Sequence`, and :class:`Scan`.
* :class:`Transducer` : represents a physical transducer device.
* :class:`Sequence` : represents pulse sequences with element-wise delays and weights.
* :class:`Scan` : represents an imaging or simulation region.
* :class:`Waveform` : represents arbitrary signals as a function of time. 
* :class:`Scatterers` : represents point scatterers in a homogeneous medium. 
* :class:`Medium` : represents arbitrary media, defined on a grid.

Simulation methods require an :class:`UltrasoundSystem` and either a :class:`Medium` or :class:`Scatterers` and produce :class:`ChannelData`.
Beamforming and post-processing methods then combine an :class:`UltrasoundSystem` and :class:`ChannelData` to produce images.

.. currentmodule:: src

Fundamental
-------------------
.. autoclass:: src.ChannelData
.. autoclass:: src.UltrasoundSystem
.. autoclass:: src.Transducer
.. autoclass:: src.Sequence
.. autoclass:: src.Scan

Simulation
------------------
.. autoclass:: src.Waveform
.. autoclass:: src.Medium
.. autoclass:: src.Scatterers

Transducers
------------
.. autoclass:: src.TransducerArray
.. autoclass:: src.TransducerConvex
.. autoclass:: src.TransducerMatrix
.. autoclass:: src.TransducerGeneric

Sequences
----------
.. autoclass:: src.SequenceRadial
.. autoclass:: src.SequenceGeneric

Scans
----------
.. autoclass:: src.ScanCartesian
.. autoclass:: src.ScanPolar
.. autoclass:: src.ScanSpherical
.. autoclass:: src.ScanGeneric
