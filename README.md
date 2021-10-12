# PVwindow: Photovoltaic Window Simulation Software

## Description

PVwindow integrates existing open source simulation codes to create a novel tool for modeling photovoltaic windows. Windows are modeled as a stack with solar irradiation impinging the outer surface with an arbitrary angle of incidence. The amount of absorption is determined in each layer. Photovoltaic layers convert some of the absorbed solar power into electrical power according to an adapted Shockley–Queisser model with a set internal quantum efficiency.[^1] The code outputs:
 - the power conversion efficiency of the photovoltaic window
 - the solar heat gain coefficient
 - The visible light transmission
 - the color of transmitted light and apparent window color window

## Implementation

PVwindow is written entirely in Python. The intensity of light that is absorbed in each layer is determined by a solution to Maxwell's equations for a stack of layers. The solution is obtained using the transfer matrix method (TMM) with an implementation that accounts for nanometer-scale layers—where interference effects are important—as well as macroscopic layers where interference effects are negligible.[^2] 

## Installation

PVwindow is a pure python code. It can be cloned into your local directory and run as long as a local installation of python is available to run it. It does depend on a few libraries, all readily available.

### Dependencies
 - numpy (https://pypi.org/project/numpy/)
 - tmm (https://pypi.org/project/tmm/)
 - colorpy (pypi repo out-dated... full python 3 version included in this repo)
 - vegas (https://pypi.org/project/vegas/)
 - pandas (https://pypi.org/project/pandas/)
 - scipy (https://pypi.org/project/scipy/)
 - matplotlib (https://pypi.org/project/matplotlib/)

## Usage

The code is organized such that a simulation script can be constructed use two types of objects: layer and stack. Layer objects store all information necessary for it to be fully described for TMM solution, its: location in the stack, thickness, and complex refractive index. A stack object defines the order of layers in the set of layers defining the window. Properties of the photovoltaic window can be processed once the stack is constructed, including:
- visible light transmission (VLT)
- solar heat gain coefficient (SHGC)
- power conversion efficiencey (PCE)
- cell temperature

### Example




[^1]: https://pubs.acs.org/doi/full/10.1021/acsenergylett.9b01316
[^2]: Code available here: https://github.com/sbyrnes321/tmm. Theory described here: https://arxiv.org/abs/1603.02720
