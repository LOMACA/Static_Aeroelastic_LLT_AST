#################################################################################

# LLT Results Display Module

################################################################################

import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from llt import LiftingLineTheory
from tkinter import filedialog

class LLTDisplay(ttk.Frame):

    """
    This class takes care of displaying the outputs of the aerodynamic calculations
    """

    def __init__(self, parent, isa_display, airfoil_aerodynamic_data_display, *args, **kwargs):

        super().__init__(parent, *args, **kwargs)
        self.isa_display = isa_display
        self.airfoil_display = airfoil_aerodynamic_data_display

        self.create_input_widgets()
        self.bind_input_clear_events()

        ttk.Button(self, text="Compute", command=self.compute_values).grid(row=10, columnspan=2, padx=5, pady=5)

        self.create_result_display()
        self.create_plot_frame()
        self.results = None

    def create_input_widgets(self):

        """
        This function creates the input widgets for the user to input
        the necessary parameters for the LLT to successfully run
        """

        self.v_inf_entry = ttk.Entry(self)
        self.AoA_entry = ttk.Entry(self)
        self.c_root_entry = ttk.Entry(self)
        self.c_tip_entry = ttk.Entry(self)
        self.b_entry = ttk.Entry(self)
        self.n_entry = ttk.Entry(self)
        self.a0_entry = ttk.Entry(self)
        self.alpha0_entry = ttk.Entry(self)

        labels = [
            ("Free Stream Velocity (m/s):", self.v_inf_entry),
            ("Angle of Attack (single value or range with delimiter -) (°):", self.AoA_entry),
            ("Root Chord Length (m):", self.c_root_entry),
            ("Tip Chord Length (m):", self.c_tip_entry),
            ("Wing Span (m):", self.b_entry),
            ("Number of Spanwise Sections:", self.n_entry),
            ("Airfoil lift curve slope (rad):", self.a0_entry),
            ("Zero lift angle of attack (°):", self.alpha0_entry),
        ]

        for row, (label_text, entry_widget) in enumerate(labels):
            ttk.Label(self, text=label_text).grid(row=row, column=0, padx=5, pady=2, sticky='w')
            entry_widget.grid(row=row, column=1, padx=5, pady=2)

    def export_all_plots(self):

        """
        This function exports all plots to the selected folder if the export option is enabled
        """

        if not self.results:
            messagebox.showwarning("Warning", "Please compute the results first.")
            return

        if not self.export_folder:
            messagebox.showerror("Export Error", "No export folder selected. Please select a folder.")
            return

        for key in ["Lift", "Bending Moment", "Downwash", "Circulation", "Induced Drag", "Lift Coefficient"]:
            self.update_plot(key)

    def bind_input_clear_events(self):

        """
        This function binds the clear_parameter function to a change in input parameter by the user
        """

        input_widgets = [
            self.v_inf_entry,
            self.AoA_entry,
            self.c_root_entry,
            self.c_tip_entry,
            self.b_entry,
            self.n_entry,
            self.a0_entry,
            self.alpha0_entry,
        ]

        for widget in input_widgets:
            widget.bind("<KeyRelease>", self.clear_results)


    def clear_results(self, event=None):

        """
        This function clears the output field
        """

        self.result_display_text.delete("1.0", tk.END)

    def create_result_display(self):

        """
       This function creates a field to display the aerodynamic results in
       """

        self.result_display_frame = ttk.LabelFrame(self, text="Results")
        self.result_display_frame.grid(row=12, column=0, columnspan=2, padx=5, pady=5)

        self.result_display_text = tk.Text(self.result_display_frame, wrap=tk.WORD, width=50, height=15)
        self.result_display_text.grid(row=13, column=0, columnspan=2, padx=5, pady=5)

    def create_plot_frame(self):

        """
        This function creates a plot frame to plot the aerodynamic results in
        """

        self.plot_frame = ttk.Frame(self)
        self.plot_frame.grid(row=0, column=2, rowspan=14, padx=5, pady=5, sticky="nsew")

        self.canvas = FigureCanvasTkAgg(plt.figure(), master=self.plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        button_frame = ttk.Frame(self)
        button_frame.grid(row=0, column=3, rowspan=14, padx=5, pady=5, sticky="nsew")

        ttk.Button(button_frame, text="Lift Distribution", command=lambda: self.update_plot("Lift")).pack(fill=tk.X)
        ttk.Button(button_frame, text="Moment Distribution", command=lambda: self.update_plot("Bending Moment")).pack(fill=tk.X)
        ttk.Button(button_frame, text="Downwash Distribution", command=lambda: self.update_plot("Downwash")).pack(fill=tk.X)
        ttk.Button(button_frame, text="Circulation Distribution", command=lambda: self.update_plot("Circulation")).pack(fill=tk.X)
        ttk.Button(button_frame, text="Induced Drag Distribution", command=lambda: self.update_plot("Induced Drag")).pack(fill=tk.X)
        ttk.Button(button_frame, text="Lift Coefficient Distribution", command=lambda: self.update_plot("Lift Coefficient")).pack(
            fill=tk.X)

        self.export_option = tk.BooleanVar()
        self.export_checkbox = ttk.Checkbutton(button_frame, text="Export Plots", variable=self.export_option)
        self.export_checkbox.pack(fill=tk.X)

        self.export_folder = None
        self.select_folder_button = ttk.Button(button_frame, text="Select Export Folder",
                                               command=self.select_export_folder)

        ttk.Button(button_frame, text="Export All Plots", command=self.export_all_plots).pack(fill=tk.X)

        self.select_folder_button.pack(fill=tk.X)
        self.select_folder_button.config(state=tk.NORMAL)

    def select_export_folder(self):

        """
        This function opens a dialog to select the folder where exported plots will be saved
        """

        folder = filedialog.askdirectory(title="Select Folder for Exported Plots")
        if folder:
            self.export_folder = folder
            messagebox.showinfo("Folder Selected", f"Export folder: {folder}")

    def compute_values(self):

        """
        This function calls the LLT calculation function and computes the aerodynamic
        results from the user inputs
        """

        try:
            self.v_inf = float(self.v_inf_entry.get())
            self.c_root = float(self.c_root_entry.get())
            self.c_tip = float(self.c_tip_entry.get())
            self.b = float(self.b_entry.get())
            self.n = int(self.n_entry.get())
            self.a0 = float(self.a0_entry.get())
            self.alpha0 = float(self.alpha0_entry.get()) * np.pi / 180

            alpha = self.AoA_entry.get().strip()

            if '-' in alpha and alpha.count('-') == 1 and alpha[0] != '-':
                start, end = map(float, alpha.split('-'))
                alpha_values = list(range(int(start), int(end) + 1))
            else:
                alpha_values = [float(alpha)]

            print("Computed AoA values:", alpha_values)

        except ValueError:
            messagebox.showerror("Error", "Invalid input. Please enter valid numeric values.")
            return

        rho = self.isa_display.rho
        tip_cl0 = self.airfoil_display.tip_cl0
        tip_cd0 = self.airfoil_display.tip_cd0
        root_cl0 = self.airfoil_display.root_cl0
        root_cd0 = self.airfoil_display.root_cd0

        lifting_line = LiftingLineTheory(self.v_inf, self.c_root, self.c_tip, self.b, self.n, root_cl0, root_cd0,
                                         self.a0, self.alpha0, rho)

        self.results = lifting_line.compute(alpha_values)

        results_text = ""

        for result in self.results:
            results_text += f"Angle of Attack: {result['alpha']}\n"
            results_text += f"CL: {result['Total CL']}\n"
            results_text += f"CD: {result['CD']}\n"
            results_text += f"Lift: {result['Total Lift']:.2f} N\n"
            results_text += f"Drag: {result['Drag']:.2f} N\n"
            results_text += f"Wing Lift Slope: {result['aw']:.2f}\n\n"

            self.aw = result['aw']
            self.L = result['Total Lift']
            self.y = result['y']
            self.AR = result['AR']
            self.TR = result['TR']
            self.S = result['S']
            self.v_inf = result['v_inf']
            self.c_root = result['c_root']
            self.c_tip = result['c_tip']
            self.b = result['b']
            self.n = result['n']
            self.L_dist = result['L_dist']
            self.L_dist_N = result['Lift']
            self.M_dist = result['Bending Moment']
            self.CL = result['Total CL']
            self.CD = result['CD']
            self.D = result['Induced Drag']
            self.cl_dist = result['Lift Coefficient']

        self.result_display_text.insert(tk.END, results_text)

    def update_plot(self, result_key):

        """
        This function updates the plot window with the aerodynamic parameter chosen by the user
        """

        if not self.results:
            messagebox.showwarning("Warning", "Please compute the results first.")
            return

        fig, ax = plt.subplots(figsize=(5, 4))

        units = {
            "Lift": "N",  # Lift distribution
            "Bending Moment": "Nm",  # Bending moment distribution
            "Downwash": "rad",  # Downwash distribution
            "Circulation": "m²/s",  # Circulation distribution
            "Induced Drag": "N/m",  # Induced drag distribution
            "Lift Coefficient": ""  # Lift coefficient distribution
        }

        for result in self.results:
            ax.plot(result['y'], result[result_key], label='LLT')

        ax.set_xlabel("y (m)")
        ax.set_ylabel(f"{result_key} ({units.get(result_key, '')})")
        #ax.set_title(f"Spanwise {result_key.replace('_', ' ').title()} Distribution")
        ax.grid(True)
        ax.legend()

        self.canvas.figure = fig
        self.canvas.draw()

        if self.export_option.get():
            if self.export_folder:
                filepath = f"{self.export_folder}/{result_key}_distribution.png"
                fig.savefig(filepath)
                messagebox.showinfo("Export Successful", f"Plot exported to {filepath}")
            else:
                messagebox.showerror("Export Error", "No export folder selected. Please select a folder.")

    def get_input_values(self):

        """
        This function retrieves the aerodynamic input values
        """

        return {
            'v_inf': self.v_inf_entry.get(),
            'AoA': self.AoA_entry.get(),
            'c_root': self.c_root_entry.get(),
            'c_tip': self.c_tip_entry.get(),
            'b': self.b_entry.get(),
            'n': self.n_entry.get(),
            'a0': self.a0_entry.get(),
            'alpha0': self.alpha0_entry.get()
        }

    def set_input_values(self, values):

        """
        This function sets the aerodynamic input values
        """

        self.v_inf_entry.delete(0, tk.END)
        self.v_inf_entry.insert(0, values.get('v_inf', self.v_inf_entry.get()))

        self.AoA_entry.delete(0, tk.END)
        self.AoA_entry.insert(0, values.get('AoA', self.AoA_entry.get()))

        self.c_root_entry.delete(0, tk.END)
        self.c_root_entry.insert(0, values.get('c_root', self.c_root_entry.get()))

        self.c_tip_entry.delete(0, tk.END)
        self.c_tip_entry.insert(0, values.get('c_tip', self.c_tip_entry.get()))

        self.b_entry.delete(0, tk.END)
        self.b_entry.insert(0, values.get('b', self.b_entry.get()))

        self.n_entry.delete(0, tk.END)
        self.n_entry.insert(0, values.get('n', self.n_entry.get()))

        self.a0_entry.delete(0, tk.END)
        self.a0_entry.insert(0, values.get('a0', self.a0_entry.get()))

        self.alpha0_entry.delete(0, tk.END)
        self.alpha0_entry.insert(0, values.get('alpha0', self.alpha0_entry.get()))

    def get_results(self):

        """
        This function collects the aerodynamic results and stores them in a dictionary
        """

        return {
            "Parameter": [
                "Lift", "Drag", "Wing Lift Slope (aw)", "Aspect Ratio (AR)", "Taper Ratio (TR)",
                "Wing Area (S)", "Free Stream Velocity (v_inf)", "Root Chord Length (c_root)",
                "Tip Chord Length (c_tip)", "Wing Span (b)", "Lift Coefficient (CL)", "Drag Coefficient (CD), Lift Coefficient Distribution"
            ],
            "Value": [
                self.L, self.D, self.aw, self.AR, self.TR, self.S, self.v_inf, self.c_root, self.c_tip, self.b, self.CL,
                self.CD, self.cl_dist
            ]
        }