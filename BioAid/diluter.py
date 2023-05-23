# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

def calculate_dilution(count_per_square: int|float, final_plate: int|float) -> None:
    """
    Calculates the dilution factor and dilution amounts for a given cell count and final plate volume.

    Args:
        count_per_square (float): The number of cells counted per square on a hemocytometer.
        final_plate (float): The desired number of cells per plate.

    Returns:
        None

    Examples:
        >>> calculate_dilution(50, 100)
        
        You should dilute by using 0.25 of cell culture with 0.75uL ddH2O 
        and then with three steps of 100uL of diluted culture and 900uL ddH2O.
        Plate with 100uL of diluted solution per plate
    """
    init_cells_ml = (count_per_square/4.0)*10**6
    final_cells_ml = final_plate*10

    dilution_factor = init_cells_ml/final_cells_ml
    dilution_factor = dilution_factor/1000

    fourth_dilution_amount = 1000/dilution_factor

    print(f"\nYou should dilute by using {fourth_dilution_amount} of cell culture with {1000-fourth_dilution_amount}uL ddH2O \nand then with three steps of 100uL of diluted culture and 900uL ddH2O.\nPlate with 100uL of diluted solution per plate")

def dilute() -> None:
    """
    Prompts the user for input and calculates the dilution factor and dilution amounts using the calculate_dilution function.

    Args:
        None

    Returns:
        None

    Raises:
        None

    Examples:
        >>> dilute()

        Welcome to the Diluter 1.0. Enjoy not doing the math yourself!

        Initial number of cells refers to the average number of cells on the 1/16th of hematocrit slide.
        For best results take the average of at least three of these squares.

        Type initial number of cells in media or enter 'x' to exit: 50

        Type in desired number of colonies per plate (100uL of media): 100

        You should dilute by using 0.25 of cell culture with 0.75uL ddH2O 
        and then with three steps of 100uL of diluted culture and 900uL ddH2O.
        Plate with 100uL of diluted solution per plate
    """

    welcome_message = "\nWelcome to the Diluter 1.0. Enjoy not doing the math yourself!\n"
    print(welcome_message)

    print("Initial number of cells refers to the average number of cells on the 1/16th of hematocrit slide.\nFor best results take the average of at least three of these squares.")

    looping = True

    while looping:
        initial_cells = input("\nType initial number of cells in media or enter 'x' to exit: ")
        if initial_cells == "x":
            looping = False
            continue
        initial_cells = float(initial_cells)
        final_cells = input("\nType in desired number of colonies per plate (100uL of media): ")
        final_cells = float(final_cells)

        calculate_dilution(initial_cells,final_cells)