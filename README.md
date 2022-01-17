# Convert cartesian to polar detector

With this tool, you can precisely convert a cartesian detector to a polar detector. 
**CartToPolarDetector** does not rely on any approximation methods like the nearest-neighbor method or interpolation. 
This tool simply overlays a cartesian detector with a polar detector grid and computes the polar pixel values 
depending on the cross-section areas with the cartesian pixels and their corresponding pixel values. 

## Installation

Run the following line:

'''
pip install -e .
'''

If you have problems with the installation you may also simply import *./src/CartToPolarDetector.py* and use it.

## How to use it

A few examples with illustrations and explanations are presented in a notebook: *./example/example_convert.ipynb*.
This short tutorial additionally requires matplotlib.

## Citation

If *ConvertToPolarDetector* is integral to a scientific publication, please cite it. A paper that describes the method that is used in this tool is submitted and will soon be found on my ORCID page: https://orcid.org/0000-0002-3639-2435

