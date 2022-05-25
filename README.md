# bio-aid
Package containing tools for genetic and genomic analysis in python in the Malkova Lab.

Includes a collection of functions from my other repositories for easier development of new projects

## Installation
The latest release of the bio-aid package can be installed through pip (https://pypi.org/project/bio-aid/) by using the following command:
```bash
pip install --upgrade bio-aid
```
bio-aid currently consists of a base module and three sub modules: `deepSeqInsH`, `MMBSearchTK`, `varaintTK`.

To easily import and access the functions inside of the bio-aid package, include the following in your `python` scripts:

```python
import BioAid as ba
```
## Diluter
This simple program calculates volumes for serial dilutions of yeast cultures, that can subsequently be used for colony plating. To access diluter, in `python` type:
```python
ba.dilute()
```
and follow instructions on the prompt.

## PopDub
Utillity for finding population doublings for the Telomere project.

To get results you will need to adjust initial and final population densities:
- set `initial_population` to the cell count on the beginning of your experiment used as inoculum.
- set `final_population_ml` to the /ml cell count at the time you want to measure population doublings.
- set `final_culture_volume` to the volume of your final culture.

Finally, run the function:
```python
calculatePopulationDoublings(initial_population, final_population_ml, final_culture_volume)
```
## RepeatSearch
This tool allows for search of imperfect repeats (Inverted and Direct) in a DNA sequence. Diagram explaining the parameters can be found in the stand-alone repo [here](https://github.com/tvarovski/RepeatSearchTools) 

```python
# Set RepeatSearch parameters
sequence="ACGT"           # This is your nucleotide sequence to be searched for repeats
inverted=True             # Sets search to Inverted (True) vs Direct (False) repeats
min_query_length=5        # Sets min length of a query sequence
max_query_length=30       # Sets max length of a query sequence
min_spacer=0              # Set min distance between query and the repeat
window_size=804           # Sets window size within which the search is confined
imperfect_homology=True   # Set True/False, to search for imperfect/perfect homologies.
min_homology=0.8          # Sets minimum homology treshold (a fraction) when imperfect_homology=True,  
fixed_errors=1            # Sets maximum number of errors (del/sub) when imperfect_homology=True (set to False or to an integer)

# To run the Search execute the following:
results_dictionary = ba.searchSequenceForRepeats(sequence=sequence,
                         min_query_length=min_query_length,
                         max_query_length=max_query_length,
                         min_spacer=min_spacer,
                         window_size=window_size,
                         imperfect_homology=imperfect_homology,
                         min_homology=min_homology,
                         fixed_errors=fixed_errors,
                         inverted=inverted
                         )

```


