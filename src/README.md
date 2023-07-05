# Class Structure
QUPS utilizes flexible abstract objects and definitions in order to preserve code re-usability. The following base classes are used to provide the majority of the functionality:

| Class | Main Properties | Main Methods |
| ------ | ------ | ------ | 
| `Transducer` | numel  | positions(), orientations() |
| `Sequence` | c0, numPulse | delays(), apodization() |
| `Scan` | size | positions(), getImagingGrid() |
| `ChannelData` | data, t0, fs | sample() |
| `Medium` | c0, rho0 | props() |
| `Scatterers` | c0, pos, amp |  |
| `Waveform` | time, samples | sample(), conv() |


A synthesis class `UltrasoundSystem` is provided to combine the `Transducer`, `Scan` and `Sequence` classes into a single definition. This is used to simulate `ChannelData` from a `Scatterer` or `Medium`, or to beamform `ChannelData` into a b-mode image. 

All of these classes provide an overloaded `plot` and/or `imagesc` method for easy visualization. 

# Data Format
QUPS objects, classes, and functions adhere to some conventions to make prototyping easier.

## Units
QUPS is unitless in time and space. Default objects are in seconds and meters, but one can easily use microseconds and millimeters for easier and quicker plotting and numerical stability. All angles are in degrees.

| Measure | Default Unit | Scalable |
| ------ | ------ | ------ |
| Position | Meter | yes |
| Time | Second | yes |
| Sound Speed | Meter/Second | yes |
| Density | Kilogram / Meter^D | yes |
| Angle | Degrees | no |

## Dimensions
All position properties and methods are defined with 3D Cartesian coordinates stored or returned in the first dimension. This includes e.g. the `Scan.positions()`, `Transducer.positions()`, and `Sequence.focus` when using a virtual source model.

Multi-dimensional arrays are described by representative characters in the "order" property of the corresponding object. For example, the data in `ChannelData.data` is described by `ChannelData.order` and the b-mode image `UltrasoundSystem.DAS` is described by `UltrasoundSystem.scan.order`.
 
| Property | Default Dimension | Swapable | 
| ------ | ------ | ------ |
| Scan.positions() | {x,y,z} in dimension 1 | no |
| Transducer.positions() | {x,y,z} in dimension 1 | no |
| ChannelData.data | Time x Receive x Transmit x F1 x F2 x ... | yes |

## Time `t = 0`

| Sequence Type | Definition | 
| ------ | ------ |
| Full-Synthetic Aperture (FSA) | Peak is centered on the firing element |
| Plane Wave (PW) | Wavefront centered on the origin (0,0,0) |
| Virtual Source (VS) | Peak is centered on the focus |

Note that this transmit sequence definition can result in data with a varying start times across transmits and a negative start time for some transmit sequence. This is explicitly supported! Simply create an array `t0` of size 1 x 1 x M, where M is the number of transmits, and use it to define the `ChannelData` e.g. with `chd = ChannelData('t0', t0, 'data', data ...)`.

# Broadcasting Support
Utilities are provided to assist in writing performant, readable code. The `sub` and `sel` utilities conveniently generate indexing expression in arbitrary dimensions. The `swapdim` utility can replace the built-in `shiftdim` and `permute` functions, and is generally more readable and more performant.

In QUPS, dimensions are used to implicitly broadcast operations, allowing you to limit memory and computational demand for simple workflows. For example, the apodization argument for the `DAS` method takes in an argument of size I1 x I2 x I3 x N x M where {I1,I2,I3} is the size of the scan, N is the number of receivers, and M is the number of transmits. This can be a huge array, easily over 100GB! However, in many cases we may only need 2 or 3 of these dimensions. 

For example, to evaluate a multi-monostatic configuration, given the FSA data, we simply create an array of size 1 x 1 x 1 x N x M and apply identity matrix weights. These weights ensure the transmitter index must match the receiver index in order to contribute to the image. No extra classes required!

```
us = UltrasoundSystem(); % create or load an UltrasoundSystem object with an FSA sequence
chd = greens(us, Scatterers('pos',[0 0 30e-3]')); % create or load a ChannelData object with FSA data
N = chd.N; % number of receivers
M = chd.M; % number of transmitters
apod = swapdim( (1:N)' == (1:M), [1,2], [4,5]); % swaps dimensions: 1 <-> 4 , 2 <-> 5
b = DAS(us, chd, 'apod', apod); % beamform using this apodization
imagesc(us.scan, b); % display the image
```
