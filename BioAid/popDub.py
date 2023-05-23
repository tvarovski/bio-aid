# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

def calculatePopulationDoublings(initial_population: int|float, final_population_ml: int|float, final_culture_volume: int|float) -> None:
    """
    Calculates the number of population doublings required to reach a final population size.

    Args:
        initial_population (int|float): The initial population size.
        final_population_ml (int|float): The final population size in milliliters.
        final_culture_volume (int|float): The final culture volume in milliliters.

    Returns:
        None: This function does not return anything, but prints the number of population doublings required to reach the final population size.
    """
    final_population = final_population_ml * final_culture_volume
    doubling_cycle_size = final_population
    doublings = 0
    while doubling_cycle_size > initial_population:
        if doubling_cycle_size < initial_population*2:
            doublings = doublings + (doubling_cycle_size-initial_population)/initial_population
            break
        doubling_cycle_size = doubling_cycle_size/2.0
        doublings += 1
        print(f"Starting from the population size of {initial_population} to the final size of {final_population} cells,\\nthere were {doublings} population doublings")