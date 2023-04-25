# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

def calculate_dilution(count_per_square, final_plate):
    init_cells_ml = (count_per_square/4.0)*10**6
    final_cells_ml = final_plate*10

    dilution_factor = init_cells_ml/final_cells_ml
    dilution_factor = dilution_factor/1000

    fourth_dilution_amount = 1000/dilution_factor

    print(f"\nYou should dilute by using {fourth_dilution_amount} of cell culture with {1000-fourth_dilution_amount}uL ddH2O \nand then with three steps of 100uL of diluted culture and 900uL ddH2O.\nPlate with 100uL of diluted solution per plate")

def dilute():

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