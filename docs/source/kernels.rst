Compute Kernels
===============
The :mod:`kern` module contains specialized compute kernels, many with CUDA or OpenCL support. 

.. currentmodule:: kern

Aperture Reduction
----------------------------
.. autofunction:: kern.cohfac
.. autofunction:: kern.dmas
.. autofunction:: kern.pcf
.. autofunction:: kern.slsc

Interpolation
----------------------------
.. autofunction:: kern.interpd
.. autofunction:: kern.rayPaths
.. autofunction:: kern.wbilerpg
.. autofunction:: kern.wbilerp
.. autofunction:: kern.wsinterpd2
.. autofunction:: kern.wsinterpd
.. autofunction:: kern.xiaolinwu_k_scaled

Signal Processing
----------------------------
.. autofunction:: kern.convd
.. autofunction:: kern.das_spec
.. autofunction:: kern.pwznxcorr

Sound Speed
----------------------------
.. autofunction:: kern.globalAverageC
.. autofunction:: kern.msfm
