# Thruster Design Utility
[![Documentation Status](https://readthedocs.org/projects/thrusterdesign/badge/?version=latest)](http://thrusterdesign.readthedocs.io/en/latest/?badge=latest)
[![Build status](https://travis-ci.org/runphilrun/TDU.svg?style=flat-square)](https://travis-ci.org/runphilrun/TDU)
[![GitHub issues](https://img.shields.io/github/issues/runphilrun/TDU.svg)](https://github.com/runphilrun/TDU/issues)
[![PyPI](https://img.shields.io/pypi/v/tdu.svg?style=flat-square)](https://pypi.python.org/pypi/tdu)
[![GitHub stars](https://img.shields.io/github/stars/runphilrun/TDU.svg)](https://github.com/runphilrun/TDU/stargazers)
[![GitHub forks](https://img.shields.io/github/forks/runphilrun/TDU.svg)](https://github.com/runphilrun/TDU/network)
[![GitHub license](https://img.shields.io/badge/license-Apache%202-blue.svg)](https://raw.githubusercontent.com/runphilrun/TDU/master/LICENSE.md)


> Tool to aid in design of small monopropellant thrusters. *[View this project on STEMN](http://stemn.com/projects/thruster-design-tool)*

## Scope
The purpose of this script is to make it easy to see how tweaks in nozzle geometry, propellants, chamber conditions, and ambient conditions affect nozzle performance in order to find the optimal design solution. This script does not perform any optimization.

## Features
* Simulate the performance of a thruster in space for a given set of parameters and output performance metrics.
* Output simulation data in a useful manner.
* Allow the user to easily tweak parameters.

## Usage
### Generating an input file
TDU loads propellant properties, inlet conditions, and nozzle geometry from a specially formatted tab-delimited text file with the extension `*.tdu`.

In general, the format of an input file is as follows:
```
PropellantNameString
Gamma    MolecularWeight
T0    P0
NumberofNodes
xLocation0    radius0
xlocation1    radius1
    .            .
    .            .
    .            .
xlocationN    radiusN
```

### Running the script
* Open `tdu.m` in Matlab 2014 or newer.
* Specify the desired input file as the value of `filein`. (For example, `filein='sample.tdu';`)
* To show additional plots and print verbose actions and data to the command line, set `debug=true;`.
* Run `tdu`
Mach number, temperature ratio, and pressure ratio at the exit of the nozzle agree with [NASA Report 1135](http://www.nasa.gov/sites/default/files/734673main_Equations-Tables-Charts-CompressibleFlow-Report-1135.pdf) for air at 1 atm with an area ratio of 2.005.


> *If you encounter any bugs, please report them in the [Issue Tracker](https://github.com/runphilrun/ThrusterDesign/issues)!*
