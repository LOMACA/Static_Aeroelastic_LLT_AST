######################################################################

# ISA Display Module for visualization in the GUI

#####################################################################


import tkinter as tk
from tkinter import ttk
from ISA import ISAComputation

class ISADisplay(ttk.Frame):

    """
    This class takes care of visualizing the ISA computations in the GUI
    """

    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        ttk.Label(self, text="Enter Altitude (km):").grid(row=0, column=0, padx=5, pady=5)
        self.altitude_entry = ttk.Entry(self)
        self.altitude_entry.grid(row=0, column=1, padx=5, pady=5)

        compute_button = ttk.Button(self, text="Compute", command=self.compute_isa)
        compute_button.grid(row=1, columnspan=2, padx=5, pady=5)

        self.result_display_text = tk.Text(self, wrap=tk.WORD, width=30, height=10)
        self.result_display_text.grid(row=2, column=0, columnspan=2, padx=5, pady=5)

    def compute_isa(self):

        """
        This function calls the ISA calculation function and displays the results in the GUI tab
        """

        try:
            altitude = float(self.altitude_entry.get())
        except ValueError:
            tk.messagebox.showerror("Error", "Please enter a valid number for altitude.")
            return

        isa_computation = ISAComputation(altitude)
        results = isa_computation.compute()

        self.rho = results['density']
        self.p = results['pressure']

        self.result_display_text.delete("1.0", tk.END)
        self.result_display_text.insert(tk.END, f"Temperature: {results['temperature']:.2f} K\n")
        self.result_display_text.insert(tk.END, f"Pressure: {results['pressure']:.2f} Pa\n")
        self.result_display_text.insert(tk.END, f"Density: {results['density']:.4f} kg/m³\n")

    def get_input_values(self):

        """
        This function collects the user input value of the ISA module (flight altitude)
        """

        return {
            'altitude': self.altitude_entry.get()
        }

    def set_input_values(self, values):

        """
        This function sets the user input value of the ISA module (flight altitude)
        """

        self.altitude_entry.delete(0, tk.END)
        self.altitude_entry.insert(0, values.get('altitude', self.altitude_entry.get()))
