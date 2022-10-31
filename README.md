# QUPS: Quick Ultrasound Processing &amp; Simulation
## Description
QUPS (pronounced "CUPS") is an abstract, lightweight, readable tool for prototyping pulse-echo ultrasound systems and algorithms. It provides a flexible, high-level representation of transducers, pulse sequences, imaging regions, and scattering media. It provides accelerated implementations of common signal processing functions for pulse-echo ultrasound data as well as cutting edge beamforming algorithms. QUPS can interface with multiple other Ultrasound simulation and processing tools including [k-Wave](http://www.k-wave.org/index.php), [MUST](https://www.biomecardio.com/MUST/documentation.html), [FieldII](https://www.field-ii.dk/) and [USTB](https://www.ustb.no/).

## Features

- Flexible:
    - 2D and 3D implementation
    - Arbitrary transducer positions and orientations
    - Arbitrary transmit waveform, delays, and apodization
    - Arbitrary pixel locations and beamforming apodization
    
- Performant:
    - Hardware acceleartion via GPU or multi-threaded CPU
    - Memory efficient data types
    - Beamform a 1024 x 1024 image for 256 x 256 transmits/receives in < 2 seconds (RTX 3070)
    - Batch simulations locally via [`parcluster`](https://www.mathworks.com/help/parallel-computing/parcluster.html) and scale to a cluster via the [parallel server](https://www.mathworks.com/help/matlab-parallel-server/) toolbox

- Modular:
    - Transducer, pulse sequence, pulse waveform, scan region etc. each defined separately
    - Easily compare waveform to waveform or transducer to transducer
    - Switch between [MUST](https://www.biomecardio.com/MUST/documentation.html), [FieldII](https://www.field-ii.dk/), or [k-Wave](http://www.k-wave.org/index.php) without redefining parameters
    - Export or import data from [USTB](https://www.ustb.no/) or [Verasonics](https://verasonics.com/vantage-advantage/) structures

- Intuitive:
    - Native MATLAB semantics with tab auto-completion
    - Overloaded `plot` and `imagesc` functions for data inspection
    - Documentation via `help` and `doc`

## Quick Start
1. Start MATLAB R2020b or later and run the setup script
```
>> cd qups; % enter the project folder
>> setup parallel CUDA; % setup the environment with any available acceleration
```

2. Create an ultrasound system and point scatterer to simulate
```
>> xdc = TransducerArray.L11_5v(); % simulate a Verasonics L11-5v transducer
>> seq = Sequence('type', 'FSA', 'numPulse', xdc.numel); % full synthetic-aperture pulse sequence
>> scan = ScanCartesian('xb', [-10e-3, 10e-3], 'zb', [25e-3 35e-3]); % set the image boundaries
>> us = UltrasoundSystem('xdc', xdc, 'sequence', seq, 'scan', scan); % create a system description
>> scat = Scatterers('pos', [0 0 30e-3]'); % a single point scatterer at 20mm depth
```

3. Display the setup
```
>> [us.scan.dx, us.scan.dz] = deal(us.lambda/4); % set the pixel resolution
>> figure;
>> plot(us); 
>> hold on;
>> plot(scat, 'cx');
```

![](fig/README/geometry.png)

4. Simulate and beamform the channel data
```
>> chd = greens(us, scat); % create channel data
>> b = DAS(us, chd); % beamform the data
>> bim = mod2db(b); % envelope detection and log-compression
```

5. Display the channel data 
```
>> figure;
>> h = imagesc(chd);
>> colormap jet; colorbar;
>> animate(h, mod2db(chd.data), 'loop', false);
```
![](fig/README/channel_data.gif)

6. Display the B-mode image
```
>> figure;
>> imagesc(us.scan, bim, max(bim(:)) + [-60 0]); % plot the image with 60dB dynamic range
>> colormap gray; colorbar;
>> title('B-mode image');
```

![](fig/README/point-target.png)



## Documentation
QUPS is documented within MATLAB. To see all the available classes, use `help ./src` or `doc ./src` from within the QUPS folder. Use `help` or `doc` on any class or method with `help classname` or `help classname.methodname` e.g. `help UltrasoundSystem.DAS`. 

## Compatibility
QUPS targets MATLAB R2020b and later. While it may work for older versions of MATLAB, you may get strange errors that don't appear in later versions. QUPS does minimal error checking in order to maintain flexibility.

If you have trouble, please submit an [issue](https://github.com/thorstone25/qups/issues).

## Extensions
All extensions to QUPS are optional, but must be installed separately from their respective sources.

| Extension | Description | Installation |
| ------ | ------ | ------ |
| [FieldII](https://www.field-ii.dk/)   | point scatterer simulator | `addpath path/to/fieldII`|
| k-Wave([base](http://www.k-wave.org/index.php), [extension](http://www.k-wave.org/forum/topic/alpha-version-of-kwavearray-off-grid-sources)) | distributed medium simulator | `addpath path/to/kWave, addpath path/to/kWaveArray` |
| [MUST](https://www.biomecardio.com/MUST/documentation.html)  | point scatterer simulator | `addpath path/to/MUST` (see issues[#2](https://github.com/thorstone25/qups/issues/2))|
| CUDA([Linux](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html),[Windows](https://docs.nvidia.com/cuda/cuda-installation-guide-microsoft-windows/index.html)) | hardware acceleration| see [CUDA Extension](####CUDA-Extension) |

### CUDA Extension
For CUDA to work, `nvcc` must succesfully run from the MATLAB environment. If a Nvidia GPU is available and `setup CUDA cache` completes with no warnings, you're all set! If you have difficulty getting nvcc to work in MATLAB, you may need to figure out which environment paths are required for _your_ CUDA installation. Running `setup CUDA` will attempt to do this for you, but may fail if you have an unexpected installation.

#### Linux
If you can run `nvcc` from a terminal or command-line interface per [CUDA installation instructions](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html), then set the `CUDA_PATH` environmental variable within MATLAB by running `setenv('CUDA_PATH', YOUR_CUDA_PATH);` prior to running `setup CUDA`. You can run `which nvcc` within a terminal to locate the CUDA_PATH installation directory. For example, if `which nvcc` returns `/opt/cuda/bin/nvcc`, then run `setenv('CUDA_PATH', 'opt/cuda');`. If this procedure does not work for you, please submit an [issue](https://github.com/thorstone25/qups/issues).

#### Windows
On Windows you must set the path for both CUDA and the _correct_ MSVC compiler for C/C++. Start a PowerShell terminal within Visual Studio. Run `echo %CUDA_PATH%` to find the base CUDA_PATH and run `echo %VCToolsInstallDir%` to find the MSVC path. Then, in MATLAB, set these paths with `setenv('CUDA_PATH', YOUR_CUDA_PATH); setenv('VCToolsInstallDir', YOUR_MSVC_PATH);`. Finally, run `setup CUDA`. From here the proper paths should be added. If this procedure does not work for you, please submit an [issue](https://github.com/thorstone25/qups/issues).
