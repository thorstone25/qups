# QUPS: Quick Ultrasound Processing &amp; Simulation

QUPS (pronounced "CUPS") is intended to be an abstract, lightweight, readable tool for prototyping simulation and processing of pulse-echo ultrasound data. It provides a flexible, high-level representation of transducers, scattering media, and accelerated implementations of common signal processing functions. Currently, QUPS can use multiple simulation tools as a backend including [k-Wave](http://www.k-wave.org/index.php), [MUST](https://www.biomecardio.com/MUST/documentation.html), and [FieldII](https://www.field-ii.dk/) and can read/write to a [USTB](https://www.ustb.no/) format.

## Getting Started
The easiest way to get started is to open and run [`example.mlx`](example.mlx) or [`example_.m`](example_.m) and interact with the simulation and beamforming examples. There are plenty of comments highlighting the different methods supported by the interface. You will need to separately download the simulator packages that you wish to use. Don't forget to add them to your path!

## Documentation
QUPS is (partially) internally documented following MATLAB conventions. This means you can use `help` or `doc` on any class or method with `help classname` or `help classname/methodname` or `help classname.methodname`. To see all the available classes, from the qups folder, use `help src` or `doc src`.

## Compatibility
QUPS targets MATLAB R2020b and later. While it may work for older versions of MATLAB, you may get strange errors that don't appear in later versions. QUPS does minimal error checking for compatibility in order to maintain flexibility.

If you have trouble, please submit an [issue](https://github.com/thorstone25/qups/issues).

### Extensions
All extensions to QUPS are optional, but must be installed separately from their respective sources.

| Extension | Installation |
| ------ | ------ |
| [FieldII](https://www.field-ii.dk/)   | add folder to path with `addpath`|
| [k-Wave](http://www.k-wave.org/index.php) | add folders for kWave and [kWaveArray](http://www.k-wave.org/forum/topic/alpha-version-of-kwavearray-off-grid-sources) to path with `addpath` |
| [MUST](https://www.biomecardio.com/MUST/documentation.html)  | add folder to path with `addpath` (see issues[#2](https://github.com/thorstone25/qups/issues/2))|
| CUDA([Linux](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html),[Windows](https://docs.nvidia.com/cuda/cuda-installation-guide-microsoft-windows/index.html)) | `!nvcc -V` must run from MATLAB |

#### CUDA Extension
If you can run `setup cache` in MATLAB, you're all set! If you have difficulty getting nvcc to work in MATLAB, you may need to figure out which environment paths are required for _your_ CUDA installation. Running `setup CUDA` will attempt to do this for you, but may fail if you have an unexpected installation.

##### Linux
If you can run `nvcc` from a terminal or command-line interface per the CUDA installation instructions, try launching MATLAB from the same [terminal](https://www.mathworks.com/help/matlab/ref/matlablinux.html). Use `getenv PATH` to display the PATH environmental variable. Try running `setup cache` to check if nvcc is working. Try the same procedure, but launching MATLAB the way you normally do. If you note any discrepancies between the paths, you can try adding each path NEW_PATH from the working setup to the system path for your normal setup in MATLAB via `setenv('PATH', fullfile(getenv('PATH'), pathsep, NEW_PATH))`. If this procedure does not work for you, please submit an issue.

##### Windows
On Windows you must include the path for both CUDA and the _correct_ MSVC compiler for C/C++. Start a PowerShell terminal within Visual Studio. Run `echo %CUDA_PATH%` to find the base CUDA_PATH and run `echo %VCToolsInstallDir%` to find the MSVC path.
Then, in MATLAB, set these paths with `setenv('CUDA_PATH', YOUR_CUDA_PATH); setenv('VCToolsInstallDir', YOUR_MSVC_PATH);`. Then run `setup CUDA`. From here the proper paths should be added. If they are not, please submit an issue.
