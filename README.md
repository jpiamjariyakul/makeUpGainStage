# A novel make-up gain stage for the software-based Moog 4-pole audio filter

Author: JAYPAT PIAMJARIYAKUL

Part of the University of Bristol's EENGM8889 "Individual Research Project 4" unit

## Prerequisites
MATLAB R2016 (or higher) must be installed.

To generate the VST audio plugin, the Audio Toolbox and Visual Studio C++ compiler 2017 (or higher) must be installed.

## Simulate Moog filter and make-up gain
In `main.m`, specify the required sampling frequency `Fs`, and the filter cutoff frequency `Fc` and feedback gain `gainFeedbk`. No further inputs are required for the make-up gain stage. Change the type of input audio signal by changing `flag` and its associated values. Run the file. Titles and captions in plots must be edited manually in-code.

Note: all the included `.m` files must be in the same directory as `main.m`.

## Generate VST plugin
In the MATLAB console, run `generateAudioPlugin vst_fullSystem.m`. The VST file should show in the working directory. Run VST with a DAW.

To generate a Win32 plugin (for compatibility), run `generateAudioPlugin --win32 vst_fullSystem.m`.