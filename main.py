import tkinter as tk
from tkinter import messagebox
import pint
import numpy as np

def dynamic_viscosity_air(T):
    return (2.791 * 10**-7 * T**0.7355) * ureg('N*s/m^2')
    
def cunningham_coefficient(dp):
    mfp = 68 * ureg('nanometer')
    return 1 + mfp/dp * (2.514 + 0.800 * np.exp(-0.55* (dp/mfp)))

def electrical_mobility(dp, C):
    e = 1.6 * 10 **(-19) * ureg('coulomb')
    mu = dynamic_viscosity_air(25+273.15)
    return (e*C)/(3*np.pi*mu*dp)
# Create a unit registry
ureg = pint.UnitRegistry()

rho = 1.293 * ureg('kg/m^3') # density of air
lambd = 6.73 * 10**-8 * ureg('meter') #mean free path
k = 1.38 * 10**-23 *ureg('J/K') # boltzmann constant

# Function to handle the submit button click
def submit():
    try:

        valid_units = {
            'Particle Diameter': ['meter']
        }
        # Read and parse input values with units
        inputs = {
        'Particle Diameter': ureg(entry1.get()),
        }
        for key, quantity in inputs.items():
            try:
                unit = str(quantity.to_base_units().units)
                if unit not in valid_units[key]:
                    raise pint.errors.UndefinedUnitError(f"Invalid unit '{unit}' for {key}")
            except (ValueError, AttributeError, pint.errors.UndefinedUnitError) as e:
                messagebox.showerror("Invalid input",
                                      f"Please enter valid floating point numbers with units. Error: {e}")
                break


        # Convert to SI units
        dp_si = inputs['Particle Diameter'].to_base_units()
        #Calculate the cunningham coefficient
        Cc = cunningham_coefficient(dp_si)
        #Format the Cunningham coefficient
        Cc_str = f'{Cc:.2f}'
        Z_p_cm = electrical_mobility(dp_si,Cc).to('cm^2/(volt*s)')
        Z_p_m = electrical_mobility(dp_si,Cc).to('m^2/(volt*s)')

        # Show the messagebox with formatted output
        messagebox.showinfo("Calculation", f'Cunningham coefficient: {Cc_str}, Electrical mobility: {round(Z_p_cm,4)}\n or {round(Z_p_m,8)}\n for D_p = {dp_si}')

    except (ValueError, pint.errors.UndefinedUnitError) as e:
        messagebox.showerror("Invalid Input", f"Please enter valid floating point numbers with units. Error: {e}")

# Create the main window
root = tk.Tk()
root.title("Cunningham coefficient calculator")

# Create labels and entry widgets
tk.Label(root, text="Enter particle diameter to calculate the corresponding Cunningham slip correction coefficient").grid(row=0, columnspan=2)

tk.Label(root, text="Particle Diameter (e.g., 5 nm):").grid(row=1, column=0)
entry1 = tk.Entry(root)
entry1.grid(row=1, column=1)

# tk.Label(root, text="Aerosol Temperature (e.g., 300 K):").grid(row=2, column=0)
# entry2 = tk.Entry(root)
# entry2.grid(row=2, column=1)

# tk.Label(root, text="Flowrate (e.g., 10 L/min):").grid(row=3, column=0)
# entry3 = tk.Entry(root)
# entry3.grid(row=3, column=1)

# tk.Label(root, text="Tube Length (e.g., 2 m):").grid(row=4, column=0)
# entry4 = tk.Entry(root)
# entry4.grid(row=4, column=1)

# tk.Label(root, text="Tube Diameter (e.g., 7 cm):").grid(row=5, column=0)
# entry5 = tk.Entry(root)
# entry5.grid(row=5, column=1)

# tk.Label(root, text="Number 6 (e.g., 50 mm):").grid(row=6, column=0)
# entry6 = tk.Entry(root)
# entry6.grid(row=6, column=1)

# Create submit and cancel buttons
tk.Button(root, text="Calculate", command=submit).grid(row=7, column=0)
tk.Button(root, text="Exit", command=root.quit).grid(row=7, column=1)

# Run the main event loop
root.mainloop()
