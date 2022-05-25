# bio-aid
Package containing tools for genetic and genomic analysis in python in the Malkova Lab.

Includes a collection of functions from my other repositories for easier development of new projects

## Installation
The latest release of the bio-aid package can be installed through pip (https://pypi.org/project/bio-aid/) by using the following command:
```bash
pip install --upgrade bio-aid
```
## Basic Usage
bio-aid currently consists of a base module and three sub modules: `deepSeqInsH`, `MMBSearchTK`, `varaintTK`.

To easily import and access the functions inside of the bio-aid package, include the following in your `python` scripts:

```python
import BioAid as ba
```
### Diluter
This simple program calculates volumes for serial dilutions of yeast cultures, that can subsequently be used for colony plating. To access diluter, in `python` type:
```python
ba.dilute()
```
and follow instructions on the prompt.

### PopDub
Utillity for finding population doublings for the Telomere project.

To get results you will need to adjust initial and final population densities:
- set `initial_population` to the cell count on the beginning of your experiment used as inoculum.
- set `final_population_ml` to the /ml cell count at the time you want to measure population doublings.
- set `final_culture_volume` to the volume of your final culture.

Finally, run the function:
```python
calculatePopulationDoublings(initial_population, final_population_ml, final_culture_volume)
```
