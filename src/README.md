# Class Structure
QUPS utilizes flexible, abstract objects and definitions in order to preserve code re-usability. The following base classes and methods are used to provide the majority of the functionality:

| Class | Main Properties | Main Methods |
| ------ | ------ | ------ | 
| `Transducer` | numel  | positions(), orientations() |
| `Sequence` | c0, numPulse | delays(), apodization() |
| `Scan` | order, size | positions() |
| `UltrasoundSystem` | xdc, seq, scan | DAS() |
| `ChannelData` | order, data, t0, fs | sample() |
| `Medium` | c0, rho0 | Sampled(), props() |
| `Scatterers` | c0, pos, amp |  |
| `Waveform` | time, samples, fs | sample(), conv() |

All of these classes provide an overloaded `plot` and/or `imagesc` method for easy visualization. 

The synthesis class `UltrasoundSystem` is provided to combine the `Transducer`, `Scan` and `Sequence` classes into a single comprehensive object. This can then be used for simulation e.g. to simulate `ChannelData` from a `Scatterer` or `Medium`, or for processing e.g. to beamform `ChannelData` into a b-mode image.

## Custom Definition Classes (Transducer, Sequence, Scan)
The Transducer, Sequence, and Scan classes are designed to be easily overloaded to maximize customizability. 

To create a custom `Transducer`, inherit from the Transducer class and define methods for the Abstract methods. For example, the following template creates a custom transducer class `TransducerCustom`, by defining the `positions()` and `orientations()` methods.
```
classdef TransducerCustom < Transducer
    properties
        ... class properties
    end
    methods
        function p = positions(xdc), ...
            ... compute element positions (3 x N)
        end
        function [az, el, normal] = orientations(xdc), ...
            ... compute element orientations (1 x N, 1 x N, 3 x N)
        end
        ... other class methods
    end
end
```

To create a custom `Sequence`, simply inherit from Sequence class and override the `delays()` and `apodization()` methods, which accept a `Transducer` as an input.
```
classdef SequenceCustom < Sequence
    properties
        ... class properties
    end
    methods
        function tau = delays(xdc), ...
            ... compute transmit delays (#elems x #pulses)
        end
        function apd = apodization(xdc), ...
            ... compute transmit element weights (#elems x #pulses)
        end
        ... other class methods
    end
end
```

To create a custom `Scan`, override the `getImagingGrid()` method, and optionally rename the variables and the character encoding of the "order" property.
```
classdef ScanCustom < Scan
    properties
        order = 'ABC';
        ... other class properties
    end
    methods
        function [X, Y, Z] = getImagingGrid(scan), ...
            ... compute pixel positions in cartesian coordinates (# A-axis, # B-axis, # C-axis)
        end
        ... other class methods
    end
end
```

## Custom Processing Methods
To define custom methods that require the full system description, you can either inherit the `UltrasoundSystem` class as above, or expand the class to include your own custom methods. To expand the `UltrasoundSystem` class, create a folder called `+UltrasoundSystem` adjacent to the `UltrasoundSystem.m` classdef file, and define a function whos first input is an `UltrasoundSystem`.
```
mkdir(fullfile(currentProject().RootFolder, "src", "+UltrasoundSystem"                    ));
edit( fullfile(currentProject().RootFolder, "src", "+UltrasoundSystem", "myCustomMethod.m")); 
```
#### myCustomMethod.m
```
function [out1, out2] = myCustomMethod(us, inp1, inp2, kwargs)
% myCustomMethod - Short description for the 'help' & 'doc'
%
% [out1, out2] = myCustomMethod(us) is the first syntax example.
%
% [...] = myCustomMethod(us, inp1, 'key1', val) is the second syntax example.
%
% See also: UltrasoundSystem.myOtherMethod UltrasoundSystem.thatOtherMethod
arguments
    us UltrasoundSystem
    inp1
    inp2 = inp1 * 2
    kwargs.key1 = "optional_default"
    kwargs.key2
end

... compute out1 and out2

end
```



# Data Format
## Units
Static constructors (e.g. `TransducerArray.L12_3v()` or `Transducer.Verasonics()`) and default objects are given in SI units with the exception of angles, which are always in degrees for readability. However, QUPS is otherwise unitless - one can easily scale metrics of time/frequency, space, and mass via the `scale` method to use e.g. microseconds/megahertz and millimeters for clarity and numerical stability.

| Measure     | Default Unit       | Scalable |
| ------      | ------             | ------   |
| Position    | Meter              | yes      |
| Time        | Second             | yes      |
| Sound Speed | Meter / Second     | yes      |
| Density     | Kilogram / Meter^3 | yes      |
| Attenuation | dB / Meter / Hz    | yes      |
| Angle       | Degrees            | no       |

**Warning**: interfacing with MUST _requires_ operating in SI units due to hard-coded constants.

## Dimensions
All position properties and methods are defined with 3D Cartesian coordinates stored or returned in the first dimension. This includes e.g. the `Scan.positions()`, `Transducer.positions()`, and `Sequence.focus` (when using a virtual source model).

Multi-dimensional arrays are described by representative characters in the "order" property of the corresponding object. For example, the data property `ChannelData.data` is described by `ChannelData.order` and the b-mode image `b = UltrasoundSystem.DAS(...)` is described by `UltrasoundSystem.scan.order`.
 
| Property               | Default Dimension                         | Swapable | 
| ---------------------- | ----------------------------------------- | -------- |
| Scan.positions()       | {x,y,z} in dimension 1                    | no       |
| Transducer.positions() | {x,y,z} in dimension 1                    | no       |
| Sequence.focus         | {x,y,z} in dimension 1                    | no       | 
| ChannelData.data       | Time x Receive x Transmit x F1 x F2 x ... | yes      |

## Time `t = 0`

| Sequence Type | Definition | 
| ----------------------------- | ---------------------------------------- |
| Full-Synthetic Aperture (FSA) | Peak is centered on the firing element   |
| Plane Wave (PW)               | Wavefront centered on the origin (0,0,0) |
| Focused (FC)                  | Peak is centered on the focus            |
| Diverging (DV)                | Peak is centered on the virtual focus    |

Note that this transmit sequence definition can result in data with a varying start times across transmits, and a negative start time for some transmit sequence. This is explicitly supported! Simply create an array `t0` of length M, where M is the number of transmits, set the size to 1 x 1 x M (assuming the transmits are in the 3rd dimension of the `data`), and define the `ChannelData` with `chd = ChannelData('t0', t0, 'data', data, 'order', 'TNM', ...)`.

# Broadcasting Support
Utilities are provided to assist in writing performant, readable code. The `sub` and `sel` utilities conveniently generate indexing expressions in arbitrary dimensions. The `swapdim` utility can replace the built-in `shiftdim` and `permute` functions, and is generally more readable and more performant.

In QUPS, dimensions are used to implicitly broadcast operations, allowing you to limit memory and computational demand for simple workflows. For example, the apodization argument for the `DAS` method takes in an argument of size I1 x I2 x I3 x N x M where {I1,I2,I3} is the size of the scan, N is the number of receivers, and M is the number of transmits. This can be a huge array, easily over 100GB! However, in many cases we may only need 2 or 3 of these dimensions. 

For example, to evaluate a multi-monostatic configuration, given the FSA data, we simply create an array of size 1 x 1 x 1 x N x M and apply identity matrix weights. These weights ensure the transmitter index must match the receiver index in order to contribute to the image. No extra classes required!

```
us = UltrasoundSystem(); % create or load an UltrasoundSystem object with an FSA sequence
scat = Scatterers('pos',[0 0 30e-3]'); % a single scaterrer
chd = greens(us, scat); % create or load a ChannelData object with FSA data
N = chd.N; % number of receivers
M = chd.M; % number of transmitters
apod = swapdim( (1:N)' == (1:M), [1,2], [4,5] ); % swaps dimensions: 1 <-> 4 , 2 <-> 5
b = DAS(us, chd, 'apod', apod); % beamform using this apodization
imagesc(us.scan, b); % display the image
```
