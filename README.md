# QUPS: Quick Ultrasound Processing &amp; Simulation

QUPS is intended to be an abstract, lightweight, readable tool for simulation and processing of pulse-echo ultrasound data to support prototyping. It provides a flexible, high-level representation of transducers, scattering media, and accelerated implementations of common signal processing functions. Currently, QUPS can use multiple simulation tools as a backend including [k-Wave](http://www.k-wave.org/index.php), [MUST](https://www.biomecardio.com/MUST/documentation.html), and [FieldII](https://www.field-ii.dk/) and can read/write to a [USTB](https://www.ustb.no/) format.

## Getting Started
The easiest way to get started is to open and run [`example.mlx`](example.mlx) or [`example_.m`](example_.m) and interact with the simulation and beamforming examples. There are plenty of comments highlighting the different methods supported by the interface. You will need to separately download the simulator packages that you wish to use. Don't forget to add them to your path!

### Compatibility
QUPS targets MATLAB R2020b and later on Linux. While it may work for older versions of MATLAB, you may get strange errors that don't appear in later versions. QUPS does minimal error checking for compatibility in order to maintain flexibility.

#### Windows Compilation
Compilation on Windows has been tested with the following dependencies:
* Windows 10 x64
* Matlab R2020a
* CUDA 10.2 and associated samples
* Visual Studio 2017

During initial setup, the user will be prompted to point to the appropriate paths to assist with compilation of the underlying CUDA binaries. 

If you have trouble compiling on any platform, please submit an [issue](https://github.com/thorstone25/qups/issues).

## Documentation
QUPS is (partially) internally documented following MATLAB conventions. This means you can use `doc` on any class and `help` on any class or method with `help classname/methodname` or `help classname.methodname`.

### Structure
QUPS utilizes flexible abstract objects and definitions in order to preserve code re-usability. The following base classes are used to provide the majority of the functionality:

| Class | Main Properties | Main Methods |
| ------ | ------ | ------ | 
| `Transducer` | numel  | positions() |
| `Sequence` | numPulse | delays(), apodization() |
| `Scan` | size | getImagingGrid() |
| `ChannelData` | data, t0, fs | sample() |
| `Target` | c0, rho0 | getPropertyMap() |

All of these classes provide an overloaded `plot` or `imagesc` method for easy visualization. 

A synthesis class `UltrasoundSystem` is provided to combine the `Transducer`, `Scan` and `Sequence` classes, simulate `ChannelData` from a `Target`, or beamform `ChannelData` into a b-mode image, which is just a regular MATLAB array ;) - you can use `Scan/imagesc` to display it.

### Data Format
QUPS objects, classes, and functions adhere to some conventions to make prototyping easier.

#### Units

| Measure | Unit | 
| ------ | ------ |
| Position | Meters |
| Time | Seconds |
| Angle | Degrees |

#### Dimensions

| Property | Standard | 
| ------ | ------ |
| Position | {x,y,z} in dimension 1 |
| rf-data | time x receive x transmit |

#### Time `t = 0`

| Sequence Type | Definition | 
| ------ | ------ |
| Full-Synthetic Aperture (FSA) | Peak is centered on the firing element |
| Plane Wave (PW) | Wavefront centered on the origin (0,0,0) |
| Virtual Source (VS) | Peak is centered on the focus |

Note that this transmit sequence definition is likely to result in data with varying start times across multiple transmits and a negative start time i.e. _before_ time 0. This is explicitly supported! Simply create a definition of `t0` of size 1 x 1 x M where M is the number of transmits that aligns the data.

#### Broadcasting
Utilities are provided to assist in writing readable, broadcasting code. The `sub` utility allows you to conveniently slice one dimension while the utility `swapdim`  as well as the built-in `shiftdim` and `permute` functions are useful for placing data in broadcasting dimensions.

Dimensions are used to implicitly broadcast operations, allowing you to limit memory and computational demand for simple workflows. For example, the apodization argument for the `DAS` method takes in an argument of size $I_1 \times I_2 \times I_3 \times N \times M$ where $\{I_1,I_2,I_3\}$ is the size of the scan, $N$ is the number of receivers, and $M$ is the number of transmits. This can be a huge array, easily over 100GB! However, in many cases we may only need 2 or 3 of these dimensions. 

For example, to evaluate a multi-monostatic configuration given a set of FSA data, we simply create an array of size $1 \times 1 \times 1 \times N \times M$ and apply identity matrix weights.
```
apod = shiftdim( (1:N)' == (1:M), -3); % we can leave this as a binary type!
```

