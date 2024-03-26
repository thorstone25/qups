Classes
=========
The :mod:`src` module contains all QUPS `classes`. 

* :class:`ChannelData` : represents the echo signal data and its time axis.
* :class:`Waveform` : represents arbitrary signals as a function of time. 
* :class:`Scatterers` : represents point scatterers in a homogeneous medium. 
* :class:`Medium` : represents arbitrary media, defined on a grid.
* :class:`Scan` : represents an imaging or simulation region.
* :class:`Transducer` : represents a physical transducer.
* :class:`Sequence` : represents pulse sequences.
* :class:`UltrasoundSystem` : combines the Sequece, Scan, and one or two Transducers.

Simulation methods require an :class:`UltrasoundSystem` and either a :class:`Medium` or :class:`Scatterers` and produce :class:`ChannelData`.
Beamforming and post-processing methods then combine an :class:`UltrasoundSystem` and :class:`ChannelData` to produce images.

.. currentmodule:: src    
.. autoclass:: src.Medium
.. autoclass:: src.Scatterers
.. autoclass:: src.SequenceGeneric
.. autoclass:: src.Sequence
.. autoclass:: src.SequenceRadial
.. autoclass:: src.TransducerArray
.. autoclass:: src.TransducerConvex
.. autoclass:: src.TransducerGeneric
.. autoclass:: src.Transducer
.. autoclass:: src.TransducerMatrix
.. autoclass:: src.Waveform
.. autoclass:: src.UltrasoundSystem