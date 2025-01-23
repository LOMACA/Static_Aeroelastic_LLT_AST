########################################################################################

# Structural model based on Schürmann "Konstruktiver Leichtbau"

########################################################################################

import numpy as np
import tkinter as tk
from tkinter import ttk
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import filedialog
from tkinter import messagebox

class CompositeModel(ttk.Frame):

    """
    This class takes care of the second structural model implementation
    designed for laminate structures, based on the description of Schürmann (2006)
    """

    def __init__(self, parent, isa_display, llt_display, airfoil_geometry_display, wing_geometry,
                 structural_model_torenbeek_display, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        self.isa_display = isa_display
        self.llt_display = llt_display
        self.airfoil_geometry_display = airfoil_geometry_display
        self.wing_geometry = wing_geometry
        self.structural_model_torenbeek_display = structural_model_torenbeek_display
        self.load_case = tk.StringVar()
        self.inertia_relief = tk.StringVar()
        self.material_properties = {}
        self.fiber_type = tk.StringVar(value='HT')

        self.create_input_widgets()
        self.show_load_case_specific_inputs()
        self.show_inertia_relief_inputs()
        self.create_plot_frame()

    def create_input_widgets(self):

        """
        This function creates input widgets for the user to input parameters
        necessary for the structural analysis
        """

        input_frame = ttk.LabelFrame(self, text="Inputs")
        input_frame.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")

        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=3)
        self.grid_columnconfigure(2, weight=10)
        self.grid_rowconfigure(0, weight=1)

        load_case_frame = ttk.LabelFrame(input_frame, text="Select Load Case")
        load_case_frame.grid(row=0, column=0, columnspan=2, padx=5, pady=5, sticky="ew")

        self.load_case.set("Maneuver")
        ttk.Radiobutton(load_case_frame, text="Maneuver", variable=self.load_case, value="Maneuver",
                        command=self.show_load_case_specific_inputs).grid(row=0, column=0, padx=5, pady=5)
        ttk.Radiobutton(load_case_frame, text="Gust", variable=self.load_case, value="Gust",
                        command=self.show_load_case_specific_inputs).grid(row=0, column=1, padx=5, pady=5)

        ttk.Label(input_frame, text="Design Maximum Take-Off Mass (kg):").grid(row=1, column=0, padx=2, pady=2)
        self.MTOM_entry = ttk.Entry(input_frame)
        self.MTOM_entry.grid(row=1, column=1, padx=2, pady=2)

        self.fiber_volume_fraction_skin_label = ttk.Label(input_frame, text="Fiber Volume Fraction (Skin) (%):")
        self.fiber_volume_fraction_skin_label.grid(row=2, column=0, padx=2, pady=2)
        self.fiber_volume_fraction_skin_entry = ttk.Entry(input_frame)
        self.fiber_volume_fraction_skin_entry.grid(row=2, column=1, padx=2, pady=2)

        self.fiber_volume_fraction_spars_label = ttk.Label(input_frame, text="Fiber Volume Fraction (Spars) (%):")
        self.fiber_volume_fraction_spars_label.grid(row=3, column=0, padx=2, pady=2)
        self.fiber_volume_fraction_spars_entry = ttk.Entry(input_frame)
        self.fiber_volume_fraction_spars_entry.grid(row=3, column=1, padx=2, pady=2)

        self.fiber_volume_fraction_skin_box_label = ttk.Label(input_frame, text="Fiber Volume Fraction (Skin Beam) (%):")
        self.fiber_volume_fraction_skin_box_label.grid(row=4, column=0, padx=2, pady=2)
        self.fiber_volume_fraction_skin_box_entry = ttk.Entry(input_frame)
        self.fiber_volume_fraction_skin_box_entry.grid(row=4, column=1, padx=2, pady=2)

        self.layer_thickness_label = ttk.Label(input_frame, text="Layer Thickness (mm):")
        self.layer_thickness_label.grid(row=5, column=0, padx=2, pady=2)
        self.layer_thickness_entry = ttk.Entry(input_frame)
        self.layer_thickness_entry.grid(row=5, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Select Fiber Type:").grid(row=6, column=0, padx=2, pady=2)
        fiber_type_dropdown = ttk.Combobox(input_frame, textvariable=self.fiber_type, values=["HT", "HM"],
                                           state="readonly")
        fiber_type_dropdown.grid(row=6, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Section Mass Axis Position (% chord):").grid(row=7, column=0, padx=2, pady=2)
        self.cg_position_entry = ttk.Entry(input_frame)
        self.cg_position_entry.grid(row=7, column=1, columnspan=2, padx=5, pady=5, sticky="ew")

        self.load_case_specific_frame = ttk.Frame(input_frame)
        self.load_case_specific_frame.grid(row=8, column=0, columnspan=2, padx=5, pady=5, sticky="ew")

        inertia_relief_frame = ttk.LabelFrame(input_frame, text="Inertia Relief")
        inertia_relief_frame.grid(row=10, column=0, columnspan=2, padx=5, pady=5, sticky="ew")

        self.inertia_relief_frame = ttk.Frame(input_frame)
        self.inertia_relief_frame.grid(row=11, column=0, columnspan=2, padx=5, pady=5, sticky="ew")

        self.inertia_relief.set("No")
        ttk.Radiobutton(inertia_relief_frame, text="No", variable=self.inertia_relief, value="No",
                        command=self.show_inertia_relief_inputs).grid(row=0, column=0, padx=5, pady=5)
        ttk.Radiobutton(inertia_relief_frame, text="Yes", variable=self.inertia_relief, value="Yes",
                        command=self.show_inertia_relief_inputs).grid(row=0, column=1, padx=5, pady=5)

        self.result_frame = ttk.LabelFrame(self, text="Results")
        self.result_frame.grid(row=0, column=1, padx=5, pady=5, sticky="nsew")

        self.result_scrollbar = ttk.Scrollbar(self.result_frame)
        self.result_scrollbar.grid(row=0, column=1, sticky="ns")

        self.result_display_text = tk.Text(self.result_frame, wrap=tk.WORD, width=70, height=50)
        self.result_display_text.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")

        self.result_display_text.config(yscrollcommand=self.result_scrollbar.set)
        self.result_scrollbar.config(command=self.result_display_text.yview)

        ttk.Button(self, text="Compute", command=lambda: self.compute(lift = self.llt_display.L_dist_N)).grid(
            row=10, column=0, padx=5, pady=5, sticky="ew")

    def show_inertia_relief_inputs(self):

        """
        This function creates input field to add inertia weights along the wing semi-span
        """

        for widget in self.inertia_relief_frame.winfo_children():
            widget.destroy()

        if self.inertia_relief.get() == "Yes":
            self.inertia_relief_entries = []

            ttk.Label(self.inertia_relief_frame, text="Add Inertia Relief Points:").grid(row=0, column=0, columnspan=2, padx=5, pady=5)

            def add_inertia_relief():
                current_row = len(self.inertia_relief_entries) + 1

                ttk.Label(self.inertia_relief_frame, text=f"Position {current_row} (section index):").grid(row=current_row * 2, column=0, padx=5, pady=5)
                position_entry = ttk.Spinbox(self.inertia_relief_frame, from_=0, to=len(self.llt_display.y) - 1, increment=1)
                position_entry.grid(row=current_row * 2, column=1, padx=5, pady=5)

                ttk.Label(self.inertia_relief_frame, text=f"Weight {current_row} (kg):").grid(row=current_row * 2 + 1, column=0, padx=5, pady=5)
                weight_entry = ttk.Entry(self.inertia_relief_frame)
                weight_entry.grid(row=current_row * 2 + 1, column=1, padx=5, pady=5)

                self.inertia_relief_entries.append((position_entry, weight_entry))

            ttk.Button(self.inertia_relief_frame, text="Add Inertia Relief Point", command=add_inertia_relief).grid(
                row=1, column=0, columnspan=4, padx=5, pady=10, sticky="ew")

        else:

            self.inertia_relief_entries = []

    def create_plot_frame(self):

        """
        This function creates a plot window to plot the results of the structural analysis
        """

        self.plot_frame = ttk.LabelFrame(self, text="Plots")
        self.plot_frame.grid(row=0, column=2, padx=5, pady=5, sticky="nsew")

        button_frame = ttk.Frame(self.plot_frame)
        button_frame.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")

        self.plot_frame.grid_rowconfigure(1, weight=1)
        self.plot_frame.grid_columnconfigure(0, weight=1)

        deflection_button = ttk.Button(button_frame, text="Plot Deflection", command=self.plot_deflection)
        deflection_button.grid(row=0, column=0, padx=5, pady=5)

        twist_button = ttk.Button(button_frame, text="Plot Twist", command=self.plot_twist)
        twist_button.grid(row=0, column=1, padx=5, pady=5)

        bending_stress_button = ttk.Button(button_frame, text="Plot Bending Stress", command=self.plot_bending_stress)
        bending_stress_button.grid(row=1, column=0, padx=5, pady=5)

        shear_stress_button = ttk.Button(button_frame, text="Plot Shear Stress", command=self.plot_shear_stress)
        shear_stress_button.grid(row=1, column=1, padx=5, pady=5)

        torsion_stress_button = ttk.Button(button_frame, text="Plot Torsion Stress", command=self.plot_torsion_stress)
        torsion_stress_button.grid(row=2, column=0, padx=5, pady=5)

        self.export_option = tk.BooleanVar()
        self.export_checkbox = ttk.Checkbutton(button_frame, text="Export Plots", variable=self.export_option)
        self.export_checkbox.grid(row=3, column=0, columnspan=2, padx=5, pady=5)

        self.export_folder = None
        self.select_folder_button = ttk.Button(button_frame, text="Select Export Folder",
                                               command=self.select_export_folder)
        self.select_folder_button.grid(row=4, column=0, columnspan=2, padx=5, pady=5)

        self.plot_canvas = None

    def select_export_folder(self):

        """
        This functions opens a dialog to select the folder where exported plots will be saved
        """

        folder = filedialog.askdirectory(title="Select Folder for Exported Plots")
        if folder:
            self.export_folder = folder
            messagebox.showinfo("Folder Selected", f"Export folder: {folder}")

    def plot_deflection(self):

        """
        This function plot the normal deflection over the span
        """

        if self.plot_canvas:
            self.plot_canvas.get_tk_widget().destroy()

        fig, ax = plt.subplots(figsize=(8, 6))

        spanwise_positions = self.llt_display.y[::-1]
        ax.plot(spanwise_positions, self.deflection, linestyle="-", label="Analytic Structural", color='b')

        ax.set_xlabel("y (m)")
        ax.set_ylabel("z (m)")
        #ax.set_title("Spanwise Z-Displacement")

        self.plot_canvas = FigureCanvasTkAgg(fig, self.plot_frame)
        self.plot_canvas.get_tk_widget().grid(row=1, column=0, padx=5, pady=5)
        ax.grid(True)
        self.plot_canvas.draw()

        if self.export_option.get():
            if self.export_folder:
                filepath = f"{self.export_folder}/deflection_plot.png"
                fig.savefig(filepath)
                messagebox.showinfo("Export Successful", f"Plot exported to {filepath}")
            else:
                messagebox.showerror("Export Error", "No export folder selected. Please select a folder.")

    def plot_twist(self):

        """
        This function plots the twist over the span
        """

        if self.plot_canvas:
            self.plot_canvas.get_tk_widget().destroy()

        fig, ax = plt.subplots(figsize=(8, 6))

        spanwise_positions = self.llt_display.y[::-1]
        ax.plot(spanwise_positions, self.twist_deg, linestyle="-", label="Analytic Structural", color='b')

        ax.set_xlabel("y (m)")
        ax.set_ylabel("Twist Angle (deg)")
        #ax.set_title("Spanwise Twist Angle")
        ax.grid(True)

        self.plot_canvas = FigureCanvasTkAgg(fig, self.plot_frame)
        self.plot_canvas.get_tk_widget().grid(row=1, column=0, padx=5, pady=5)
        self.plot_canvas.draw()

        if self.export_option.get():
            if self.export_folder:
                filepath = f"{self.export_folder}/twist_plot.png"
                fig.savefig(filepath)
                messagebox.showinfo("Export Successful", f"Plot exported to {filepath}")
            else:
                messagebox.showerror("Export Error", "No export folder selected. Please select a folder.")

    def plot_bending_stress(self):

        """
        This function plots the bending stress over the span
        """

        if self.plot_canvas:
            self.plot_canvas.get_tk_widget().destroy()

        spanwise_position = self.llt_display.y[::-1]
        bending_stress_spar1 = self.bending_stress_spar1_list
        bending_stress_spar2 = self.bending_stress_spar2_list
        bending_stress_upper = self.bending_stress_upper_sparcap_list
        bending_stress_lower = self.bending_stress_lower_sparcap_list

        fig, ax = plt.subplots()
        ax.plot(spanwise_position, bending_stress_spar1, label="Spar Web 1 Bending Stress")
        ax.plot(spanwise_position, bending_stress_spar2, label="Spar Web 2 Bending Stress")
        ax.plot(spanwise_position, bending_stress_upper, label="Upper Sparcap Bending Stress")
        ax.plot(spanwise_position, bending_stress_lower, label="Lower Sparcap Bending Stress")

        #ax.set_title("Spanwise Bending Stress Distribution")
        ax.set_xlabel("y (m)")
        ax.set_ylabel("Bending Stress (MPa)")
        ax.grid(True)
        ax.legend()

        self.plot_canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        self.plot_canvas.get_tk_widget().grid(row=1, column=0, padx=5, pady=5, sticky="nsew")
        self.plot_canvas.draw()

        if self.export_option.get():
            if self.export_folder:
                filepath = f"{self.export_folder}/bending_stress_plot.png"
                fig.savefig(filepath)
                messagebox.showinfo("Export Successful", f"Plot exported to {filepath}")
            else:
                messagebox.showerror("Export Error", "No export folder selected. Please select a folder.")

    def plot_shear_stress(self):

        """
        This function plots the shear stress over the span
        """

        if self.plot_canvas:
            self.plot_canvas.get_tk_widget().destroy()

        spanwise_position = self.llt_display.y[::-1]
        shear_stress_spar1 = self.shear_stress_spar1_list
        shear_stress_spar2 = self.shear_stress_spar2_list
        shear_stress_upper = self.shear_stress_upper_sparcap_list
        shear_stress_lower = self.shear_stress_lower_sparcap_list

        fig, ax = plt.subplots()
        ax.plot(spanwise_position, shear_stress_spar1, label="Spar Web 1 Shear Stress")
        ax.plot(spanwise_position, shear_stress_spar2, label="Spar Web 2 Shear Stress")
        ax.plot(spanwise_position, shear_stress_upper, label="Upper Sparcap Shear Stress")
        ax.plot(spanwise_position, shear_stress_lower, label="Lower Sparcap Shear Stress")

        #ax.set_title("Spanwise Shear Stress Distribution")
        ax.set_xlabel("y (m)")
        ax.set_ylabel("Shear Stress (MPa)")
        ax.grid(True)
        ax.legend()

        self.plot_canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        self.plot_canvas.get_tk_widget().grid(row=1, column=0, padx=5, pady=5, sticky="nsew")
        self.plot_canvas.draw()

        if self.export_option.get():
            if self.export_folder:
                filepath = f"{self.export_folder}/shear_stress_plot.png"
                fig.savefig(filepath)
                messagebox.showinfo("Export Successful", f"Plot exported to {filepath}")
            else:
                messagebox.showerror("Export Error", "No export folder selected. Please select a folder.")

    def plot_torsion_stress(self):

        """
        This function plots the torsion stress over the span
        """

        if self.plot_canvas:
            self.plot_canvas.get_tk_widget().destroy()

        spanwise_position = self.llt_display.y[::-1]
        torsion_stress_nose = self.torsion_stress_nose_list
        torsion_stress_spar = self.torsion_stress_center_list
        torsion_stress_rear = self.torsion_stress_rear_list
        torsion_stress_spar1 = self.torsion_stress_spar1_list
        torsion_stress_spar2 = self.torsion_stress_spar2_list

        fig, ax = plt.subplots()
        ax.plot(spanwise_position, torsion_stress_nose, label="Nose Torsion Stress")
        ax.plot(spanwise_position, torsion_stress_spar, label="Box Torsion Stress")
        ax.plot(spanwise_position, torsion_stress_rear, label="Rear Torsion Stress")
        ax.plot(spanwise_position, torsion_stress_spar1, label="Spar Web 1 Torsion Stress")
        ax.plot(spanwise_position, torsion_stress_spar2, label="Spar Web 2 Torsion Stress")

        #ax.set_title("Spanwise Torsion Stress Distribution")
        ax.set_xlabel("y (m)")
        ax.set_ylabel("Torsion Stress (MPa)")
        ax.grid(True)
        ax.legend()

        self.plot_canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        self.plot_canvas.get_tk_widget().grid(row=1, column=0, padx=5, pady=5, sticky="nsew")
        self.plot_canvas.draw()

        if self.export_option.get():
            if self.export_folder:
                filepath = f"{self.export_folder}/torsion_stress_plot.png"
                fig.savefig(filepath)
                messagebox.showinfo("Export Successful", f"Plot exported to {filepath}")
            else:
                messagebox.showerror("Export Error", "No export folder selected. Please select a folder.")

    def get_mtom(self):

        """
        This function retrieves the MTOM of the aircraft from inputs
        """

        return float(self.MTOM_entry.get())

    def show_load_case_specific_inputs(self):

        """
        This function displays inputs necessary to compute the load factor for two load cases
        """

        for widget in self.load_case_specific_frame.winfo_children():
            widget.destroy()

        if self.load_case.get() == "Gust":
            ttk.Label(self.load_case_specific_frame, text="Derived Gust Velocity (m/s):").grid(row=0, column=0, padx=5,
                                                                                               pady=5)
            self.gust_velocity_entry = ttk.Entry(self.load_case_specific_frame)
            self.gust_velocity_entry.grid(row=0, column=1, padx=5, pady=5)

            ttk.Label(self.load_case_specific_frame, text="Equivalent Airspeed (m/s):").grid(row=1, column=0, padx=5,
                                                                                             pady=5)
            self.equivalent_airspeed_entry = ttk.Entry(self.load_case_specific_frame)
            self.equivalent_airspeed_entry.grid(row=1, column=1, padx=5, pady=5)

    def set_material_properties(self):

        """
        This function assigns material properties based on the chosen material type
        """

        fiber_type = self.fiber_type.get()

        properties = {
            "HT": {
                "E_parallel": 230000,  # N/mm^2
                "E_perpendicular": 28000,  # N/mm^2
                "G_parallel_perpendicular": 50000,  # N/mm^2
                "nu_parallel_perpendicular": 0.23,
                "rho_f": 1.74 * 1000  # kg/m^3
            },
            "HM": {
                "E_parallel": 392000,  # N/mm^2
                "E_perpendicular": 15200,  # N/mm^2
                "G_parallel_perpendicular": 28600,  # N/mm^2
                "nu_parallel_perpendicular": 0.2,
                "rho_f": 1.81 * 1000 # kg/m^3
            }
        }

        self.material_properties = properties[fiber_type]

    def compute(self, lift):

        """
        This function performs the structural analysis
        """

        self.result_display_text.delete("1.0", tk.END)
        MTOM = float(self.MTOM_entry.get())

        results_text = ""

        fiber_type = self.fiber_type.get()

        if fiber_type == "HT":
            E_parallel = 230000 * 1e6
            E_perpendicular = 28000 * 1e6
            G_parallel_perpendicular = 50000 * 1e6
            nu_parallel_perpendicular = 0.23
            rho_f = 1.74 * 1000
        elif fiber_type == "HM":
            E_parallel = 392000 * 1e6
            E_perpendicular = 15200 * 1e6
            G_parallel_perpendicular = 28600 * 1e6
            nu_parallel_perpendicular = 0.2
            rho_f = 1.81 * 1000

        else:
            self.result_display_text.insert(tk.END, "Error: No valid fiber type selected.\n")
            return

        fiber_volume_fraction_skin = float(self.fiber_volume_fraction_skin_entry.get()) / 100
        fiber_volume_fraction_spars = float(self.fiber_volume_fraction_spars_entry.get()) / 100
        fiber_volume_fraction_spar_caps = float(self.fiber_volume_fraction_skin_box_entry.get()) / 100

        self.bending_stress_spar1_list = []
        self.bending_stress_spar2_list = []
        self.bending_stress_upper_sparcap_list = []
        self.bending_stress_lower_sparcap_list = []

        self.shear_stress_spar1_list = []
        self.shear_stress_spar2_list = []
        self.shear_stress_upper_sparcap_list = []
        self.shear_stress_lower_sparcap_list = []

        self.torsion_stress_nose_list = []
        self.torsion_stress_center_list = []
        self.torsion_stress_rear_list = []
        self.torsion_stress_spar1_list = []
        self.torsion_stress_spar2_list = []

        self.total_weight_distribution = []

        for spanwise_position in range(len(self.llt_display.y)):
            result_text = f"Spanwise Position {spanwise_position}:\n"

            if self.load_case.get() == "Maneuver":
                self.load_factor = self.calculate_maneuver_load_factor(MTOM)
                result_text += f"Load Factor: {self.load_factor:.2f}\n"
                self.total_lift = np.sum(lift[spanwise_position:])
                result_text += f"Total lift up to spanwise section {spanwise_position}: {self.total_lift:.2f} N\n"

            else:
                gust_velocity = float(self.gust_velocity_entry.get())
                equivalent_airspeed = float(self.equivalent_airspeed_entry.get())
                self.load_factor = self.calculate_gust_load_factor(MTOM, gust_velocity, equivalent_airspeed)
                result_text += f"Load Factor: {self.load_factor:.2f}\n"
                self.total_lift = np.sum(lift[spanwise_position:])
                result_text += f"Total lift up to spanwise section {spanwise_position}: {self.total_lift:.2f} N\n"

            airfoil, partition1, partition2 = self.wing_geometry.wing_geometry[spanwise_position]

            self.spar1_height, self.spar2_height = self.calculate_spar_heights(airfoil, spanwise_position)
            self.section_length = np.sqrt((self.llt_display.y[-1] - self.llt_display.y[spanwise_position]) ** 2 + (
                    airfoil[-1][2] - airfoil[0][2]) ** 2)
            self.skin_area = self.calculate_skin_area(airfoil) * self.section_length
            self.spar_thickness = self.wing_geometry.get_spar_thickness(spanwise_position)

            A_soll_spar1 = self.spar_thickness
            A_soll_spar2 = self.spar_thickness

            partition1_x = self.airfoil_geometry_display.partition1 * self.wing_geometry.get_chord_length(
                spanwise_position)
            partition2_x = self.airfoil_geometry_display.partition2 * self.wing_geometry.get_chord_length(
                spanwise_position)

            self.l_gurt = partition2_x - partition1_x
            self.skin_thickness = float(self.airfoil_geometry_display.skin_thickness_entry.get())
            self.skin_thickness_box = float(self.airfoil_geometry_display.skin_thickness_box_entry.get())

            A_soll_gurt = self.skin_thickness_box

            m_f_over_L = 800*1e-6 # Roving Feinheit in Tex (g/km)

            n_spar1 = int(A_soll_spar1 * (rho_f * fiber_volume_fraction_spars) / m_f_over_L)
            n_spar2 = int(A_soll_spar2 * (rho_f * fiber_volume_fraction_spars) / m_f_over_L)
            n_gurt = int(A_soll_gurt * (rho_f * fiber_volume_fraction_skin) / m_f_over_L)

            total_thickness_spar1 = self.spar_thickness
            total_thickness_spar2 = self.spar_thickness

            layer_thickness = float(self.layer_thickness_entry.get())/1000

            num_layers_spar1 = int(total_thickness_spar1 / layer_thickness)
            num_layers_spar2 = int(total_thickness_spar2 / layer_thickness)
            num_layers_gurt = int(self.skin_thickness_box / layer_thickness)

            result_text += f"Material: {fiber_type}\n"
            #result_text += f"Spar 1 Height: {self.spar1_height:.4f} m\n"
            #result_text += f"Spar 2 Height: {self.spar2_height:.4f} m\n"
            #result_text += f"Spar Cap Length: {l_gurt:.4f} m\n"
            #result_text += f"Number of Layers Spar 1: {num_layers_spar1}\n"
            #result_text += f"Number of Layers Spar 2: {num_layers_spar2}\n"
            #result_text += f"Number of Layers Spar Caps: {num_layers_gurt}\n"
            #result_text += f"Number Rovings Spar 1: {n_spar1}\n"
            #result_text += f"Number Rovings Spar 2: {n_spar2}\n"
            #result_text += f"Number Rovings Spar Cap: {n_gurt}\n"

            layer_orientations = [0, 45, -45, 90]

            half_layers_spar1 = num_layers_spar1 // 2
            half_layers_spar2 = num_layers_spar2 // 2
            half_layers_gurt = num_layers_gurt // 2

            # Matrix Values (EP Harz für Segelflugzeuge)

            E_m = 3150 * 1e6
            v_m = 0.37
            rho_m = 1.19 * 1000

            self.E_parallel_skin = self.calculate_youngs_modulus_parallel(E_parallel, fiber_volume_fraction_skin, E_m)
            self.E_parallel_spar = self.calculate_youngs_modulus_parallel(E_parallel, fiber_volume_fraction_spars, E_m)
            self.E_parallel_spar_cap = self.calculate_youngs_modulus_parallel(E_parallel, fiber_volume_fraction_spar_caps, E_m)

            self.E_perp_skin = self.calculate_youngs_modulus_orthogonal(E_m, v_m, fiber_volume_fraction_skin, E_perpendicular)
            self.E_perp_spar = self.calculate_youngs_modulus_orthogonal(E_m, v_m, fiber_volume_fraction_spars, E_perpendicular)
            self.E_perp_spar_cap = self.calculate_youngs_modulus_orthogonal(E_m, v_m, fiber_volume_fraction_spar_caps,
                                                                        E_perpendicular)

            self.G_skin = self.calculate_shear_modulus(fiber_volume_fraction_skin, G_parallel_perpendicular)
            self.G_spar = self.calculate_shear_modulus(fiber_volume_fraction_spars,
                                                                    G_parallel_perpendicular)
            self.G_spar_cap = self.calculate_shear_modulus(fiber_volume_fraction_spar_caps,
                                                       G_parallel_perpendicular)

            self.nu_skin = self.calculate_poissons_ratio(fiber_volume_fraction_skin, nu_parallel_perpendicular, v_m)
            self.nu_spar = self.calculate_poissons_ratio(fiber_volume_fraction_spars, nu_parallel_perpendicular, v_m)
            self.nu_spar_cap = self.calculate_poissons_ratio(fiber_volume_fraction_spar_caps, nu_parallel_perpendicular, v_m)

            Q_matrices_spar1, S_matrices_spar1 = self.create_symmetric_matrices(half_layers_spar1, num_layers_spar1, layer_orientations,
                                                                           self.E_parallel_spar, self.E_perp_spar,
                                                                           self.G_spar,
                                                                           self.nu_spar)


            Q_matrices_spar2, S_matrices_spar2 = self.create_symmetric_matrices(half_layers_spar2, num_layers_spar2, layer_orientations,
                                                                           self.E_parallel_spar, self.E_perp_spar,
                                                                           self.G_spar,
                                                                           self.nu_spar)


            Q_matrices_gurt, S_matrices_gurt = self.create_symmetric_matrices(half_layers_gurt, num_layers_gurt, layer_orientations,
                                                                         self.E_parallel_spar_cap, self.E_perp_spar_cap,
                                                                         self.G_spar_cap,
                                                                         self.nu_spar_cap)

            t_k = layer_thickness
            t_spar1 = total_thickness_spar1
            t_spar2 = total_thickness_spar2
            self.t_gurt = self.skin_thickness_box

            self.A_matrix_spar1 = self.calculate_A_matrix(Q_matrices_spar1, t_k)
            self.A_matrix_spar2 = self.calculate_A_matrix(Q_matrices_spar2, t_k)
            self.A_matrix_gurt1 = self.calculate_A_matrix(Q_matrices_gurt, t_k)
            self.A_matrix_gurt2 = self.calculate_A_matrix(Q_matrices_gurt, t_k)

            self.E_x_spar1, self.E_y_spar1, self.G_xy_spar1, self.nu_xy_spar1 = self.calculate_effective_moduli(self.A_matrix_spar1, t_spar1)

            self.E_x_spar2, self.E_y_spar2, self.G_xy_spar2, self.nu_xy_spar2 = self.calculate_effective_moduli(self.A_matrix_spar2, t_spar2)

            self.E_x_gurt, self.E_y_gurt, self.G_xy_gurt, self.nu_xy_gurt = self.calculate_effective_moduli(self.A_matrix_gurt1, self.t_gurt)


            #self.E_x_spar1 = 71000e6
            #self.E_x_spar2 = 71000e6
            #self.E_x_gurt = 71000e6
            #self.E_y_spar1 = 71000e6
            #self.E_y_spar2 = 71000e6
            #self.E_y_gurt = 71000e6
            #self.G_xy_spar1 = 27000e6
            #self.G_xy_spar2 = 27000e6
            #self.G_xy_gurt = 27000e6
            #self.nu_xy_spar1 = 0.33
            #self.nu_xy_spar2 = 0.33
            #self.nu_xy_gurt = 0.33

            result_text += f"Stiffness Matrix Spar 1:\n{self.A_matrix_spar1}\n"
            result_text += f"Stiffness Matrix Spar 2:\n{self.A_matrix_spar2}\n"
            result_text += f"Stiffness Matrix Sparcap 1:\n{self.A_matrix_gurt1}\n"
            result_text += f"Stiffness Matrix Sparcap 2:\n{self.A_matrix_gurt2}\n"

            result_text += f"Laminate Young's Modulus (Spar 1): {self.E_x_spar1/1e6:.2f} N/mm^2 (Parallel)\n"
            result_text += f"Laminate Young's Modulus  (Spar 2): {self.E_x_spar2/1e6:.2f} N/mm^2 (Parallel)\n"
            result_text += f"Laminate Young's Modulus (Spar Cap): {self.E_x_gurt/1e6:.2f} N/mm^2 (Parallel)\n"
            result_text += f"Laminate Young's Modulus(Spar 1): {self.E_y_spar1 / 1e6:.2f} N/mm^2 (Perpendicular)\n"
            result_text += f"Laminate Young's Modulus (Spar 2): {self.E_y_spar2 / 1e6:.2f} N/mm^2 (Perpendicular)\n"
            result_text += f"Laminate Shear Modulus (Spar 1): {self.G_xy_spar1 / 1e6:.2f} N/mm^2\n"
            result_text += f"Laminate Shear Modulus (Spar 2): {self.G_xy_spar2 / 1e6:.2f} N/mm^2\n"
            result_text += f"Laminate Shear Modulus (Spar Cap): {self.G_xy_gurt/1e6:.2f} N/mm^2\n"
            result_text += f"Laminate Poisson's Ratio (Spar 1): {self.nu_xy_spar1:.3f}\n"
            result_text += f"Laminate Poisson's Ratio (Spar 2): {self.nu_xy_spar2:.3f}\n"
            result_text += f"Laminate Poisson's Ratio (Spar Cap): {self.nu_xy_gurt:.3f}\n"

            self.material_density_skin = self.calculate_composite_density(fiber_volume_fraction_skin, rho_f)
            self.material_density_spars = self.calculate_composite_density(fiber_volume_fraction_spars, rho_f)

            #self.material_density_skin = 2795
            #self.material_density_spars = 36/1000

            self.spar1_length = np.sqrt((self.llt_display.y[-1] - self.llt_display.y[spanwise_position]) ** 2 + (
                    airfoil[-1][2] - airfoil[0][2]) ** 2)
            self.spar2_length = self.spar1_length

            self.spar1_weight = self.spar1_length * self.spar1_height * t_spar1 * self.material_density_spars
            self.spar2_weight = self.spar2_length * self.spar2_height * t_spar2 * self.material_density_spars

            self.skin_circumference = self.calculate_circumference(airfoil)
            self.skin_weight = self.skin_circumference * float(
                self.airfoil_geometry_display.skin_thickness_entry.get()) * self.material_density_skin

            #self.wing_weight = 0.328

            inertia_relief_weight_at_section = 0

            if self.inertia_relief.get() == "Yes":
                for position_entry, weight_entry in self.inertia_relief_entries:
                    try:
                        position = int(position_entry.get())
                        weight = float(weight_entry.get())
                    except ValueError:
                        continue

                    if spanwise_position == position:
                        inertia_relief_weight_at_section += weight

            self.total_weight = (self.spar1_weight + self.spar2_weight + self.skin_weight + inertia_relief_weight_at_section) * 9.8

            self.total_weight_distribution.append(self.total_weight)

            # Torsion

            self.torsion_moment = self.calculate_torsion_moment(spanwise_position, lift)
            self.shear_flows, self.nu = self.calculate_torsion(airfoil, spanwise_position, self.torsion_moment, self.G_xy_gurt, self.G_xy_spar1)
            self.nu_deg = self.nu * (180/np.pi)
            self.n_xs1 = float(self.shear_flows[0])
            self.n_xs2 = float(self.shear_flows[1])
            self.n_xs3 = float(self.shear_flows[2])
            self.n_xs_spar_1 = self.n_xs1 - self.n_xs2
            self.n_xs_spar_2 = self.n_xs2 - self.n_xs3
            self.torsion_shear_stresses = self.calculate_shear_stress(self.shear_flows)
            self.torsion_shear_strains = self.calculate_strain(self.torsion_shear_stresses, self.E_y_spar1)
            self.skin_thickness = float(self.airfoil_geometry_display.skin_thickness_entry.get())
            self.torsion_shear_stresses_mm = [stress / 1e6 for stress in self.torsion_shear_stresses]
            self.torsion_shear_stress_spar_1 = self.n_xs_spar_1 / self.spar_thickness
            self.torsion_shear_stress_spar_2 = self.n_xs_spar_2 / self.spar_thickness
            self.torsion_shear_strain_spar_1 = self.torsion_shear_stress_spar_1 / self.E_y_spar1
            self.torsion_shear_strain_spar_2 = self.torsion_shear_stress_spar_2 / self.E_y_spar2
            self.torsion_shear_stress_spar_1_mm = self.torsion_shear_stress_spar_1 / 1e6
            self.torsion_shear_stress_spar_2_mm = self.torsion_shear_stress_spar_2 / 1e6

            self.components = [
                {'E': self.E_x_spar1, 'A': self.spar_thickness * self.spar1_height,
                 'z': self.spar1_height, 'type': 'spar', 'height': self.spar1_height,
                 'thickness': self.spar_thickness},
                {'E': self.E_x_spar2, 'A': self.spar_thickness * self.spar2_height,
                 'z': self.spar2_height, 'type': 'spar', 'height': self.spar2_height,
                 'thickness': self.spar_thickness},
                {'E': self.E_x_gurt, 'A': self.t_gurt * self.l_gurt, 'z': self.spar1_height / 2,
                 'type': 'cap'},
                {'E': self.E_x_gurt, 'A': self.t_gurt * self.l_gurt, 'z': -self.spar2_height / 2,
                 'type': 'cap'}
            ]

            # Bending

            self.bending_moment = self.calculate_bending_moment(spanwise_position, lift) #*self.load_factor)
            self.bending_strain_distribution, self.bending_stress_distribution = self.calculate_bending_stress_strain(spanwise_position, lift)
            self.bending_stress_distribution_mpa = [stress / 1e6 for stress in self.bending_stress_distribution]
            self.bending_stress_spar_1 = self.bending_stress_distribution[0]
            self.bending_stress_spar2 = self.bending_stress_distribution[1]
            self.bending_stress_upper_sparcap = self.bending_stress_distribution[2]
            self.bending_stress_lower_sparcap = -self.bending_stress_distribution[3]
            self.bending_strain_spar_1 = self.bending_strain_distribution[0]
            self.bending_strain_spar_2 = self.bending_strain_distribution[1]
            self.bending_strain_upper_sparcap = self.bending_strain_distribution[2]
            self.bending_strain_lower_sparcap = -self.bending_strain_distribution[3]
            self.bending_stress_spar_1_mm = self.bending_stress_distribution[0]/1e6
            self.bending_stress_spar2_mm = self.bending_stress_distribution[1]/1e6
            self.bending_stress_upper_sparcap_mm = self.bending_stress_distribution[2]/1e6
            self.bending_stress_lower_sparcap_mm = -self.bending_stress_distribution[3]/1e6

            # Shear

            self.components_shear = [
                {'E': self.E_y_spar1, 'A': self.spar_thickness * self.spar1_height,
                 'z': self.spar1_height, 'type': 'spar', 'height': self.spar1_height,
                 'thickness': self.spar_thickness},
                {'E': self.E_y_spar2, 'A': self.spar_thickness * self.spar2_height,
                 'z': -self.spar2_height, 'type': 'spar', 'height': self.spar2_height,
                 'thickness': self.spar_thickness},
                {'E': self.E_y_gurt, 'A': self.t_gurt * self.l_gurt, 'z': self.spar1_height / 2,
                 'type': 'cap', 'height': self.spar1_height, 'thickness': self.t_gurt},
                {'E': self.E_y_gurt, 'A': self.t_gurt * self.l_gurt, 'z': self.spar2_height / 2,
                 'type': 'cap', 'height': self.spar2_height, 'thickness': self.t_gurt}

            ]

            self.neutral_axis = self.calculate_neutral_axis()
            #result_text += f"Neutral axis: {self.neutral_axis}\n"
            I_y_shear = self.calculate_moment_of_inertia()
            self.I_y_shear = I_y_shear
            self.I_y = self.calculate_geometric_moment_of_inertia()
            #result_text += f"Geometric Moment of Inertia: {self.I_y}\n"

            self.Q_z = self.calculate_shear_force(spanwise_position, lift, self.total_weight) # *self.load_factor, self.total_weight)
            self.S_s_list = self.calculate_static_moment()
            #result_text += f"Qz: {self.Q_z}\n"
            #result_text += f"Static Moments: {self.S_s_list}\n"
            #result_text += f"Youngs Modulus: {self.E}\n"
            self.shear_flows_shear = self.calculate_shear_flow(self.Q_z, self.E_y_spar1, self.S_s_list, I_y_shear)
            self.shear_flow_shear_spar1 = self.shear_flows_shear[0]
            self.shear_flow_shear_spar2 = self.shear_flows_shear[1]
            self.shear_flow_shear_sparcap1  =self.shear_flows_shear[2]
            self.shear_flow_shear_sparcap2 = self.shear_flows_shear[3]
            self.shear_stress_shear_spar1 = self.shear_flow_shear_spar1 / self.spar_thickness
            self.shear_stress_shear_spar2 = self.shear_flow_shear_spar2 / self.spar_thickness
            self.shear_stress_shear_sparcap1 = self.shear_flow_shear_sparcap1 / self.skin_thickness_box
            self.shear_stress_shear_sparcap2 = self.shear_flow_shear_sparcap2 / self.skin_thickness_box
            self.shear_stress_shear_spar1_mm = self.shear_stress_shear_spar1 / 1e6
            self.shear_stress_shear_spar2_mm = self.shear_stress_shear_spar2 / 1e6
            self.shear_stress_shear_sparcap1_mm  = self.shear_stress_shear_sparcap1 / 1e6
            self.shear_stress_shear_sparcap2_mm = self.shear_stress_shear_sparcap2 / 1e6
            self.shear_strain_shear_sparcap1 = self.shear_stress_shear_sparcap1 / self.E_y_gurt
            self.shear_strain_shear_sparcap2 = self.shear_stress_shear_sparcap2 / self.E_y_gurt
            self.shear_strain_shear_spar1 = self.shear_stress_shear_spar1 / self.E_y_spar1
            self.shear_strain_shear_spar2 = self.shear_stress_shear_spar2 / self.E_y_spar2

            #result_text += f"Spar 1 Height: {self.spar1_height:.2f} m\n"
            #result_text += f"Spar 2 Height: {self.spar2_height:.2f} m\n"
            #result_text += f"Spar Length: {self.spar1_length:.2f} m\n"
            #result_text += f"Skin Area: {self.skin_area:.2f} m^2\n"
            #result_text += f"Skin Circumference: {self.skin_circumference:.2f} m\n"
            #result_text += f"Material Density (Skin): {self.material_density_skin:.2f} kg/m^3\n"
            #result_text += f"Material Density (Spars): {self.material_density_spars:.2f} kg/m^3\n"
            #result_text += f"Spar 1 Weight: {self.spar1_weight:.2f} N\n"
            #result_text += f"Spar 2 Weight: {self.spar2_weight:.2f} N\n"
            #result_text += f"Skin Weight: {self.skin_weight:.2f} N\n"
            #result_text += f"Wing Weight: {self.wing_weight:.2f} N\n"
            result_text += f"Total Weight: {self.total_weight_distribution[spanwise_position]:.2f} N\n"
            #result_text += f"Young's Modulus (Skin): {self.youngs_modulus_skin/1e6:.2f} N/mm^2\n"
            #result_text += f"Young's Modulus (Spars): {self.youngs_modulus_spars/1e6:.2f} N/mm^2\n"
            #result_text += f"Shear Modulus (Skin): {self.shear_modulus_skin/1e6:.2f} N/mm^2\n"
            #result_text += f"Shear Modulus (Spars): {self.shear_modulus_spars/1e6:.2f} N/mm^2\n"
            result_text += f"Torsion Moment: {self.torsion_moment:.2f} N*m\n"
            result_text += f"Bending Moment: {self.bending_moment:.2f} N*m\n"
            result_text += f"Shear Flows due to Torsion (Nose, Spar, Rear): {self.shear_flows} N/m \n"
            result_text += f"Shear Flows due to Torsion (Sparweb 1 and 2): {self.n_xs_spar_1, self.n_xs_spar_2} N/m\n"
            result_text += f"Shear Stresses Torsion (Nose, Spar, Rear): {self.torsion_shear_stresses_mm} N/mm^2 \n"
            result_text += f"Shear Stresses Torsion (Sparweb 1, Sparweb 2): {self.torsion_shear_stress_spar_1_mm, self.torsion_shear_stress_spar_2_mm} N/mm^2\n"
            result_text += f"Shear Strains Torsion (Nose, Spar, Rear): {self.torsion_shear_strains}\n"
            result_text += f"Shear Strains Torsion (Spar 1, Spar 2): {self.torsion_shear_strain_spar_1, self.torsion_shear_strain_spar_2} \n"
            result_text += f"Warping: {self.nu_deg:.6e} deg/m\n"
            #result_text += f"Bending Strain Distribution: {self.bending_strain_distribution}\n"
            #result_text += f"Bending Stress Distribution: {self.bending_stress_distribution_mpa} N/mm^2 \n"
            result_text += f"Bending Stress Spar 1: {self.bending_stress_spar_1_mm:.2f} N/mm^2\n"
            result_text += f"Bending Stress Spar 2: {self.bending_stress_spar2_mm:.2f} N/mm^2\n"
            result_text += f"Bending Stress Upper Sparcap: {self.bending_stress_upper_sparcap_mm:.2f} N/mm^2\n"
            result_text += f"Bending Stress Lower Sparcap: {self.bending_stress_lower_sparcap_mm:.2f} N/mm^2\n"
            result_text += f"Bending Strain Spar 1: {self.bending_strain_spar_1:.6f}\n"
            result_text += f"Bending Strain Spar 2: {self.bending_strain_spar_2:.6f}\n"
            result_text += f"Bending Strain Upper Sparcap: {self.bending_strain_upper_sparcap:.6f}\n"
            result_text += f"Bending Strain Lower Sparcap: {self.bending_strain_lower_sparcap:.6f}\n"
            result_text += f"Shear Flows due to Shear (Web 1, Web 2, Upper Sparcap, Lower Sparcap): {self.shear_flows_shear} N/m \n"
            result_text += f"Shear Stress due to Shear (Web 1, Web 2, Upper Sparcap, Lower Sparcap): {self.shear_stress_shear_spar1_mm, self.shear_stress_shear_spar2_mm, self.shear_stress_shear_sparcap1_mm, self.shear_stress_shear_sparcap2_mm} N/mm^2 \n"
            result_text += f"Shear Strain due to Shear (Web 1, Web 2, Upper Sparcap, Lower Sparcap): {self.shear_strain_shear_spar1, self.shear_strain_shear_spar2, self.shear_strain_shear_sparcap1, self.shear_strain_shear_spar2} \n"

            results_text += result_text + "\n"

            self.bending_stress_spar1_list.append(self.bending_stress_spar_1_mm)
            self.bending_stress_spar2_list.append(self.bending_stress_spar2_mm)
            self.bending_stress_upper_sparcap_list.append(self.bending_stress_upper_sparcap_mm)
            self.bending_stress_lower_sparcap_list.append(self.bending_stress_lower_sparcap_mm)

            self.shear_stress_spar1_list.append(self.shear_stress_shear_spar1_mm)
            self.shear_stress_spar2_list.append(self.shear_stress_shear_spar2_mm)
            self.shear_stress_upper_sparcap_list.append(self.shear_stress_shear_sparcap1_mm)
            self.shear_stress_lower_sparcap_list.append(self.shear_stress_shear_sparcap2_mm)

            self.torsion_stress_nose_list.append(self.torsion_shear_stresses_mm[0])
            self.torsion_stress_center_list.append(self.torsion_shear_stresses_mm[1])
            self.torsion_stress_rear_list.append(self.torsion_shear_stresses_mm[2])
            self.torsion_stress_spar1_list.append(self.torsion_shear_stress_spar_1_mm)
            self.torsion_stress_spar2_list.append(self.torsion_shear_stress_spar_2_mm)

            self.result_display_text.insert(tk.END, results_text)

        self.deflection = self.calculate_deflection(self.E_x_spar1, lift)
        self.twist = self.calculate_twist_distribution()
        self.twist_deg = [twist * (180/np.pi) for twist in self.twist]

        result_text = ""

        result_text += f"Deflection: {self.deflection} m \n"
        result_text += f"Twist: {self.twist_deg} deg \n"

        results_text += result_text + "\n"
        self.result_display_text.insert(tk.END, results_text)

    def calculate_transformed_Q_matrix(self, E_parallel, E_perpendicular, G_parallel_perpendicular, nu_parallel_perpendicular,
                                       alpha):

        """
        This function calculates the element stiffness matrix in the CLT
        """

        alpha_rad = np.radians(alpha)

        Q11 = E_parallel / (1 - nu_parallel_perpendicular ** 2)
        Q22 = E_perpendicular / (1 - nu_parallel_perpendicular ** 2)
        Q12 = nu_parallel_perpendicular * E_perpendicular / (1 - nu_parallel_perpendicular ** 2)
        Q66 = G_parallel_perpendicular

        cos_alpha = np.cos(alpha_rad)
        sin_alpha = np.sin(alpha_rad)
        cos2_alpha = cos_alpha ** 2
        sin2_alpha = sin_alpha ** 2
        cos4_alpha = cos2_alpha ** 2
        sin4_alpha = sin2_alpha ** 2
        sin2_alpha_2 = np.sin(2*alpha_rad)**2

        Q_bar11 = Q11 * cos4_alpha + Q22 * sin4_alpha + 0.5 * (Q12 + 2 * Q66) * sin2_alpha_2
        Q_bar22 = Q11 * sin4_alpha + Q22 * cos4_alpha + 0.5 * (Q12 + 2 * Q66) * sin2_alpha_2
        Q_bar12 = Q12 + 0.25 * (Q11 + Q22 - 2*Q12 - 4*Q66) * sin2_alpha_2
        Q_bar66 = Q66 + 0.25 * (Q11 + Q22 - 2*Q12 - 4*Q66) * sin2_alpha_2
        Q_bar16 = -0.5 * ((Q11 + Q22 - 2*Q12 - 4*Q66) * sin2_alpha - (Q11 - Q12 - 2*Q66)) * np.sin(2*alpha_rad)
        Q_bar26 = -0.5 * ((Q22 - Q12 - 2*Q66) - (Q11 + Q22 - 2*Q12 - 4*Q66) * sin2_alpha) * np.sin(2*alpha_rad)

        Q_bar = np.array([
            [Q_bar11, Q_bar12, Q_bar16],
            [Q_bar12, Q_bar22, Q_bar26],
            [Q_bar16, Q_bar26, Q_bar66]
        ])

        return Q_bar

    def calculate_transformed_S_matrix(self, E_parallel, E_perpendicular, G_parallel_perpendicular, nu_parallel_perpendicular,
                                       alpha):

        """
        This function calculates the element strain matrix
        """

        alpha_rad = np.radians(alpha)

        S11 = 1/E_parallel
        S22 = 1/E_perpendicular
        S66 = 1/G_parallel_perpendicular
        S12 = nu_parallel_perpendicular/E_parallel

        cos_alpha = np.cos(alpha_rad)
        sin_alpha = np.sin(alpha_rad)
        cos2_alpha = cos_alpha ** 2
        sin2_alpha = sin_alpha ** 2
        sin3_alpha = sin2_alpha * sin_alpha
        cos3_alpha = cos2_alpha * cos_alpha
        cos4_alpha = cos2_alpha ** 2
        sin4_alpha = sin2_alpha ** 2
        sin2_alpha_2 = np.sin(2*alpha_rad)**2
        cos2_alpha_2 = np.cos(2*alpha_rad)**2

        S_bar11 = cos4_alpha/E_parallel + sin4_alpha/E_perpendicular + 0.25 * (S66 - 2*S12) * sin2_alpha_2
        S_bar22 = sin4_alpha/E_parallel + cos4_alpha/E_perpendicular + 0.25 * (S66 - 2*S12) * sin2_alpha_2
        S_bar66 = cos2_alpha_2/G_parallel_perpendicular + (S11 + S22 + 2*S12) * sin2_alpha_2
        S_bar12 = 0.25 * (S11 + S22 - S66) * sin2_alpha_2 - S12 * (sin4_alpha + cos4_alpha)
        S_bar16 = -(2/E_perpendicular + 2*S12 - S66) * sin3_alpha * cos_alpha + (2/E_parallel + 2*S12 - S66) * cos3_alpha * sin_alpha
        S_bar26 = -(2/E_perpendicular + 2*S12 - S66) * cos3_alpha * sin_alpha + (2/E_parallel + 2*S12 - S66) * sin3_alpha * cos_alpha

        S_bar = np.array([
            [S_bar11, S_bar12, S_bar16],
            [S_bar12, S_bar22, S_bar26],
            [S_bar16, S_bar26, S_bar66]
        ])

        return S_bar

    def create_symmetric_matrices(self, num_layers_half, num_layers_spar1, layer_orientations, E_parallel, E_perpendicular,
                                  G_parallel_perpendicular, nu_parallel_perpendicular):

        """
        This function ensures that the input layer angles lead to symmetric stiffness matrices built
        in the CLT
        """

        Q_matrices = []
        S_matrices = []

        for i in range(num_layers_half):
            orientation = layer_orientations[i % len(layer_orientations)]
            Q_bar = self.calculate_transformed_Q_matrix(E_parallel, E_perpendicular, G_parallel_perpendicular,
                                                        nu_parallel_perpendicular, orientation)
            S_bar = self.calculate_transformed_S_matrix(E_parallel, E_perpendicular, G_parallel_perpendicular,
                                                        nu_parallel_perpendicular, orientation)
            Q_matrices.append(Q_bar)
            S_matrices.append(S_bar)


        if num_layers_half * 2 < num_layers_spar1:

            orientation = layer_orientations[0]
            Q_bar = self.calculate_transformed_Q_matrix(E_parallel, E_perpendicular, G_parallel_perpendicular,
                                                        nu_parallel_perpendicular, orientation)
            S_bar = self.calculate_transformed_S_matrix(E_parallel, E_perpendicular, G_parallel_perpendicular,
                                                        nu_parallel_perpendicular, orientation)
            Q_matrices.append(Q_bar)
            S_matrices.append(S_bar)


        for i in range(num_layers_half - 1, -1, -1):
            Q_matrices.append(Q_matrices[i])
            S_matrices.append(S_matrices[i])

        return Q_matrices, S_matrices

    def calculate_A_matrix(self, Q_matrices, layer_thickness):

        """
        This function computes the laminate stiffness matrix in the CLT
        """

        A_matrix = np.zeros((3, 3))

        for Q_bar in Q_matrices:
            A_matrix[0, 0] += (Q_bar[0, 0] * layer_thickness)
            A_matrix[0, 1] += (Q_bar[0, 1] * layer_thickness)
            A_matrix[0, 2] += (Q_bar[0, 2] * layer_thickness)
            A_matrix[1, 0] += (Q_bar[1, 0] * layer_thickness)
            A_matrix[1, 1] += (Q_bar[1, 1] * layer_thickness)
            A_matrix[1, 2] += (Q_bar[1, 2] * layer_thickness)
            A_matrix[2, 0] += (Q_bar[2, 0] * layer_thickness)
            A_matrix[2, 1] += (Q_bar[2, 1] * layer_thickness)
            A_matrix[2, 2] += (Q_bar[2, 2] * layer_thickness)

        return A_matrix

    def calculate_effective_moduli(self, A_matrix, t):

        """
        This function calculates the laminate mechanical properties
        """

        A_inv = np.linalg.inv(A_matrix)

        A_inv_11 = A_inv[0, 0]
        A_inv_22 = A_inv[1, 1]
        A_inv_66 = A_inv[2, 2]
        A_inv_12 = A_inv[0, 1]

        E_x = 1 / (A_inv_11 * t)
        E_y = 1 / (A_inv_22 * t)
        G_xy = 1 / (A_inv_66 * t)

        nu_xy = -A_inv_12 / A_inv_22
        nu_yx = -A_inv_12 / A_inv_11

        return E_x, E_y, G_xy, nu_xy

    def calculate_youngs_modulus_parallel(self, E_f_parallel, fiber_volume_fraction, E_m):

        """
        This function calculates the layer Young's Modulus parallel
        to the fiber direction based on the chosen material type
        """

        return E_f_parallel * fiber_volume_fraction + E_m * (1-fiber_volume_fraction)

    def calculate_youngs_modulus_orthogonal(self, E_m, v_m, fiber_volume_fraction, E_f_perp):

        """
        This function calculates the layer Young's Modulus perpendicular
        to the fiber direction based on the chosen material type
        """

        return (E_m / (1 - v_m**2)) * (1 + 0.85 * fiber_volume_fraction**2) / \
         ((1 - fiber_volume_fraction)**1.25 + (E_m / ((1 - v_m**2) * E_f_perp)) * fiber_volume_fraction)

    def calculate_shear_modulus(self, fiber_volume_fraction, G_f):

        """
        This function calculates the layer Shear Modulus
        """

        G_m = 1019*1e6 # N/m^2

        return G_m * (1 + 0.4 * fiber_volume_fraction ** 0.5) / (
                    (1 - fiber_volume_fraction) ** 1.45 + (G_m / G_f) * fiber_volume_fraction)

    def calculate_poissons_ratio(self, fiber_volume_fraction, v_f, v_m):

        """
        This function calculates the layer Poisson's Ratio
        """

        return fiber_volume_fraction * v_f + (1-fiber_volume_fraction) * v_m

    def calculate_maneuver_load_factor(self, MTOM):

        """
        This function calculates the maneuver load factor according to CS23 Easy Access Rules
        """

        MTOM = float(self.MTOM_entry.get())
        MTOW = MTOM * 9.8
        MTOW_lb = MTOW / 4.448
        n_man = 2.1 + (24000 / (MTOW_lb + 10000))
        return n_man

    def calculate_gust_load_factor(self, MTOM, gust_velocity, equivalent_airspeed):

        """
        This function calculates the gust load factor according to CS23 Easy Access Rules
        """

        U_de = float(self.gust_velocity_entry.get())
        V = float(self.equivalent_airspeed_entry.get())
        rho0 = 1.225
        rho = self.isa_display.rho
        S = self.llt_display.S
        g = 9.8
        W = float(self.MTOM_entry.get()) * g
        a = self.llt_display.aw
        mgc = self.wing_geometry.calculate_mgc()
        mug = (2 * (W / S)) / (rho * mgc * a * g)
        kg = (0.88 * mug) / (5.3 + mug)
        n_gust = 1 + (kg * rho0 * U_de * V * a) / (2 * (W / S))
        return n_gust


    def compute_lift(self, spanwise_position):

        """
        This function calculates the spanwise lift acting at the respective station
        """

        L_dist = self.llt_display.L_dist
        lift = np.sum(L_dist[spanwise_position:])

        return lift

    def calculate_spar_heights(self, airfoil, spanwise_position):

        """
        This function calculates the heights of the spar webs
        """

        partition1_x = self.airfoil_geometry_display.partition1 * self.wing_geometry.get_chord_length(spanwise_position)
        partition2_x = self.airfoil_geometry_display.partition2 * self.wing_geometry.get_chord_length(spanwise_position)

        partition1_z_upper = self.find_z_at_x(partition1_x, airfoil, upper=True)
        partition1_z_lower = self.find_z_at_x(partition1_x, airfoil, upper=False)
        partition2_z_upper = self.find_z_at_x(partition2_x, airfoil, upper=True)
        partition2_z_lower = self.find_z_at_x(partition2_x, airfoil, upper=False)

        spar1_height = abs(partition1_z_upper - partition1_z_lower)
        spar2_height = abs(partition2_z_upper - partition2_z_lower)

        return spar1_height, spar2_height

    def find_z_at_x(self, x_pos, coordinates, upper=True):

        """
        This function finds the corresponding z-coordinate value at the chordwise position x
        """

        if upper:
            filtered_coords = [point for point in coordinates if point[1] > np.mean([p[1] for p in coordinates])]
        else:
            filtered_coords = [point for point in coordinates if point[1] <= np.mean([p[1] for p in coordinates])]

        if not filtered_coords:
            return 0

        closest_index = np.argmin([abs(point[0] - x_pos) for point in filtered_coords])
        return filtered_coords[closest_index][2]

    def calculate_skin_area(self, coordinates):

        """
        This function calculates area covered by the wing skin
        """

        upper_surface = [coord for coord in coordinates if coord[1] > np.mean([p[1] for p in coordinates])]
        lower_surface = [coord for coord in coordinates if coord[1] <= np.mean([p[1] for p in coordinates])]

        def calculate_surface_area(surface):
            x = np.array([coord[0] for coord in surface])
            z = np.array([coord[2] for coord in surface])
            area = np.trapz(z, x)
            return area

        upper_area = calculate_surface_area(upper_surface)
        lower_area = calculate_surface_area(lower_surface)

        return abs(upper_area) + abs(lower_area)

    def calculate_circumference(self, coordinates):

        """
        This function calculates circumference of each cross-section
        """

        upper_surface = [coord for coord in coordinates if coord[1] > np.mean([p[1] for p in coordinates])]
        lower_surface = [coord for coord in coordinates if coord[1] <= np.mean([p[1] for p in coordinates])]

        def calculate_surface_circumference(surface):
            x = np.array([coord[0] for coord in surface])
            z = np.array([coord[2] for coord in surface])
            circumference = np.sum(np.sqrt(np.diff(x) ** 2 + np.diff(z) ** 2))
            return circumference

        upper_circumference = calculate_surface_circumference(upper_surface)
        lower_circumference = calculate_surface_circumference(lower_surface)

        return abs(upper_circumference) + abs(lower_circumference)

    def calculate_composite_density(self, fiber_volume_fraction, rho_f):

        """
        This function calculates the density of the composite material from
        fiber and matrix densities
        """

        rho_m = 1190
        return rho_f * fiber_volume_fraction + rho_m * (1 - fiber_volume_fraction)

    def calculate_torsion(self, airfoil, spanwise_position, torsion_moment, G_skin, G_spars):

        """
        This function implements the torsion calculation for a three-cell cross-section
        described by Schürmann
        """

        G_skin = self.G_xy_gurt
        G_spars = self.G_xy_spar1

        nose_thickness = float(self.airfoil_geometry_display.skin_thickness_entry.get())
        spar_thickness = float(self.airfoil_geometry_display.spar_thickness_entry.get())

        spar1_height, spar2_height = self.calculate_spar_heights(airfoil, spanwise_position)

        spar_length = self.wing_geometry.get_chord_length(spanwise_position) * (
                    self.airfoil_geometry_display.partition2 - self.airfoil_geometry_display.partition1)
        nose_length = self.calculate_circumference(airfoil) * self.airfoil_geometry_display.partition1
        trailing_edge_length = self.calculate_circumference(airfoil) * (1 - self.airfoil_geometry_display.partition2)

        A_m1 = self.calculate_cell_area(airfoil, 1, spanwise_position)
        A_m2 = self.calculate_cell_area(airfoil, 2, spanwise_position)
        A_m3 = self.calculate_cell_area(airfoil, 3, spanwise_position)

        def equations(shear_flows):
            n_xs1, n_xs2, n_xs3 = shear_flows

            eq1 = torsion_moment - (2 * A_m1 * n_xs1 + 2 * A_m2 * n_xs2 + 2 * A_m3 * n_xs3)

            g_prime_1 = (n_xs1 / (2 * A_m1 * G_skin) * nose_length / nose_thickness) + (n_xs1 - n_xs2) / (
                        2 * A_m1 * G_spars * spar1_height / spar_thickness)
            g_prime_2 = (n_xs2 / (2 * A_m2 * G_skin) * spar_length / nose_thickness) + (n_xs2 - n_xs1) / (
                        2 * A_m2 * G_spars * spar1_height / spar_thickness) + (n_xs2 - n_xs3) / (
                                    2 * A_m2 * G_spars * spar2_height / spar_thickness)
            g_prime_3 = (n_xs3 / (2 * A_m3 * G_skin) * trailing_edge_length / nose_thickness) + (n_xs2 - n_xs3) / (
                        2 * A_m3 * G_spars * spar2_height / spar_thickness)

            eq2 = g_prime_1 - g_prime_2
            eq3 = g_prime_2 - g_prime_3

            return [eq1, eq2, eq3]

        initial_guess = np.array([1.0, 1.0, 1.0])

        shear_flows = fsolve(equations, initial_guess)

        n_xs1 = shear_flows[0]

        ds_over_t_sum = (nose_length / nose_thickness) + (spar1_height / spar_thickness) + (
                    spar2_height / spar_thickness) + (trailing_edge_length / nose_thickness)

        nu = n_xs1 / (2 * A_m1 * G_skin) * ds_over_t_sum

        return shear_flows, nu

    def calculate_cell_area(self, airfoil, cell_id,spanwise_position):

        """
        This function calculates the area of each cell out of the three cells constituting the
        cross-section
        """

        if cell_id == 1:

            spar1_height, _ = self.calculate_spar_heights(airfoil, spanwise_position)
            nose_length = self.calculate_circumference(airfoil) * self.airfoil_geometry_display.partition1
            return 0.5 * nose_length * spar1_height
        elif cell_id == 2:
            spar1_height, spar2_height = self.calculate_spar_heights(airfoil, spanwise_position)
            spar_length = self.wing_geometry.get_chord_length(spanwise_position) * (
                        self.airfoil_geometry_display.partition2 - self.airfoil_geometry_display.partition1)
            return spar_length * (spar1_height + spar2_height) / 2
        elif cell_id == 3:
            _, spar2_height = self.calculate_spar_heights(airfoil, spanwise_position)
            trailing_edge_length = self.calculate_circumference(airfoil) * (
                        1 - self.airfoil_geometry_display.partition2)
            return 0.5 * trailing_edge_length * spar2_height
        else:
            return 0

    def calculate_torsion_moment(self, spanwise_position, lift):

        """
        This function calculates the torsional moment acting at each cross-section
        """

        total_lift = np.sum(lift[spanwise_position:])
        distance_to_cg = float(self.cg_position_entry.get()) / 100 * self.wing_geometry.get_chord_length(spanwise_position)
        torsion_moment = total_lift * distance_to_cg

        return torsion_moment

    def calculate_spar_width(self, spanwise_position):

        """
        This function calculates the distance between the two spar webs
        """

        partition1_x = self.airfoil_geometry_display.partition1 * self.wing_geometry.get_chord_length(spanwise_position)
        partition2_x = self.airfoil_geometry_display.partition2 * self.wing_geometry.get_chord_length(spanwise_position)
        spar_width = partition2_x - partition1_x
        return spar_width

    def calculate_neutral_axis(self):

        """
        This function determines the position of the neutral axis
        """

        numerator = sum(E_i * A_i * z_i for component in self.components for E_i, A_i, z_i in [(component['E'], component['A'], component['z'])])
        denominator = sum(E_i * A_i for component in self.components for E_i, A_i in [(component['E'], component['A'])])
        z_neutral = numerator / denominator
        return z_neutral

    def calculate_moment_of_inertia(self):

        """
        This function calculates the second moments of inertia of the spar components
        """

        z_neutral = self.calculate_neutral_axis()
        I_y = 0

        for component in self.components:
            E_i = component['E']
            A_i = component['A']
            z_i = component['z']

            if component['type'] == 'spar':
                h = component['height']
                t = component['thickness']
                I_y_spar = (t * h ** 3) / 12
                I_y += E_i * I_y_spar

            elif component['type'] == 'cap':
                I_y_cap = A_i * (z_i - z_neutral) ** 2
                I_y += E_i * I_y_cap

        return I_y

    def calculate_weighted_E_eff(self):

        """
        This function calculates the effective Young's Modulus of the spar
        """

        E_weighted_sum = 0
        I_sum = 0

        for component in self.components:
            E_i = component['E']
            I_i = self.calculate_geometric_moment_of_inertia_component(component)
            E_weighted_sum += E_i * I_i
            I_sum += I_i

        if I_sum == 0:
            return 0
        E_eff = E_weighted_sum / I_sum

        return E_eff

    def calculate_geometric_moment_of_inertia_component(self, component):

        """
        This function calculates the second moments of inertia of the spar components
        """

        I_y_component  = []

        if component['type'] == 'spar':
            h = component['height']
            t = component['thickness']

            I_y_component = (t * h ** 3) / 12

        elif component['type'] == 'cap':
            A_cap = component['A']
            z_i = component['z']
            z_neutral = self.calculate_neutral_axis()

            I_y_component = A_cap * (z_i - z_neutral) ** 2

        return I_y_component

    def calculate_deflection(self, E, lift):

        """
        This function calculates the bending deflection
        """

        deflection = [0] * len(self.llt_display.y)
        slope = [0] * len(self.llt_display.y)

        print(
            f"{'Position':>10} | {'Bending Moment':>16} | {'Moment of Inertia (I_y)':>26} | {'Slope':>10} | {'Deflection':>12}")

        for i in range(1, len(self.llt_display.y)):

            I_y = self.calculate_geometric_moment_of_inertia()

            #lift = self.compute_total_lift(i, self.load_factor)
            #lift = np.sum(lift_dist[i:])*self.load_factor

            M_b = self.calculate_bending_moment(i, lift)

            dx = self.llt_display.y[i] - self.llt_display.y[i - 1]

            slope[i] = slope[i - 1] + (M_b / (E * I_y)) * dx

            deflection[i] = deflection[i - 1] + slope[i] * dx

            print(
                f"{self.llt_display.y[i]:>10.2f} m | {M_b:>16.2f} Nm | {I_y:>26.6e} m^4 | {slope[i]:>10.6f} | {deflection[i]:>12.6f} m")

        return deflection

    def calculate_geometric_moment_of_inertia(self):

        """
        This function calculates the second moment of inertia
        of the entire spar
        """

        t_spar = self.spar_thickness
        h_spar1 = self.spar1_height
        h_spar2 = self.spar2_height

        A_spar_cap = self.skin_thickness_box * self.l_gurt
        z_spar_cap = h_spar1 / 2

        I_spar1 = (t_spar * h_spar1 ** 3) / 12
        I_spar2 = (t_spar * h_spar2 ** 3) / 12

        I_gurt = 2 * A_spar_cap * z_spar_cap ** 2

        I_y = I_spar1 + I_spar2 + I_gurt

        return I_y

    def calculate_twist_distribution(self):

        """
        This function calculates the wing twist angle from the twist rate
        """

        y = self.llt_display.y[::-1]

        twist = [0] * len(y)

        for i in range(1, len(y)):

            dx = y[i] - y[i - 1]

            twist[i] = twist[i - 1] + self.nu * dx

        return twist

    def calculate_bending_stress_strain(self, spanwise_position, lift):

        """
        This function calculates the bending stress
        """

        z_neutral = self.calculate_neutral_axis()

        I_y = self.calculate_moment_of_inertia()

        M_b = self.calculate_bending_moment(spanwise_position, lift) #*self.load_factor)

        strain_distribution = []
        stress_distribution = []

        for component in self.components:
            E_i = component['E']
            z_i = component['z']

            strain = (M_b / I_y) * (z_i - z_neutral)
            strain_distribution.append(strain)

            stress = E_i * strain
            stress_distribution.append(stress)

        return strain_distribution, stress_distribution

    def calculate_bending_moment(self, spanwise_position, lift):

        """
        This function calculates the bending moment at each spanwise station
        """

        y_position = self.llt_display.y[spanwise_position]
        total_lift = np.sum(lift[spanwise_position:])
        M = (total_lift - self.total_weight) * y_position
        return M

    def calculate_strain(self, stresses, E):

        """
        This function calculates strains from stresses through Hooke's Law
        """

        strains = [stress / E for stress in stresses]
        return strains

    def calculate_shear_stress(self, torsion_flows):

        """
        This function calculates the shear stress
        """

        skin_thickness = float(self.airfoil_geometry_display.skin_thickness_entry.get())
        stresses = [torsion_flow / skin_thickness for torsion_flow in torsion_flows]
        return stresses

    def calculate_shear_force(self, spanwise_position, lift, total_weight):

        """
        This function calculates the shear force at each spanwise station
        """

        total_lift = np.sum(lift[spanwise_position:])
        Q_z = total_lift - total_weight

        return Q_z

    def calculate_static_moment(self):

        """
        This function calculates the static moment at each spanwise station
        """

        z_neutral = self.calculate_neutral_axis()
        S_s_list = []

        for component in self.components_shear:
            A_i = component['A']
            z_i = component['z']

            if component['type'] == 'spar':
                h = component['height']
                t = component['thickness']
                A_spar = h * t
                S_s = A_spar * (z_i - z_neutral)
                S_s_list.append(S_s)

            elif component['type'] == 'cap':
                S_s = A_i * (z_i - z_neutral)
                S_s_list.append(S_s)

        return S_s_list

    def calculate_shear_flow(self, Q_z, E,  S_s_list, I_y_shear):

        """
        This function calculates the shear flows at each spanwise station
        """

        shear_flows = []
        for S_s in S_s_list:
            n_xs = (-Q_z * E * S_s) / I_y_shear
            shear_flows.append(n_xs)
        return shear_flows

    def get_input_values(self):

        """
        This function collects the input parameters
        """

        values = {
            'load_case': self.load_case.get(),
            'MTOM': self.MTOM_entry.get(),
            'inertia_relief': self.inertia_relief.get(),
            'cg_position': self.cg_position_entry.get(),
            'fiber_volume_fraction_skin': self.fiber_volume_fraction_skin_entry.get(),
            'fiber_volume_fraction_spars': self.fiber_volume_fraction_spars_entry.get(),
            'fiber_volume_fraction_beam': self.fiber_volume_fraction_skin_box_entry.get(),
            'layer_thickness': self.layer_thickness_entry.get()
        }

        if self.load_case.get() == "Gust":
            values.update({
                'gust_velocity': self.gust_velocity_entry.get(),
                'equivalent_airspeed': self.equivalent_airspeed_entry.get()
            })

        if self.inertia_relief.get() == "Yes" and hasattr(self, 'inertia_relief_entries'):
            inertia_positions = []
            inertia_weights = []

            for position_entry, weight_entry in self.inertia_relief_entries:
                try:
                    position = float(position_entry.get())
                    weight = float(weight_entry.get())
                    inertia_positions.append(position)
                    inertia_weights.append(weight)
                except ValueError:

                    continue

            values.update({
                'inertia_relief_positions': inertia_positions,
                'inertia_relief_weight': inertia_weights
            })

        return values

    def set_input_values(self, values):

        """
        This function sets the saved input values
        """

        self.load_case.set(values.get('load_case', 'Maneuver'))
        self.MTOM_entry.delete(0, tk.END)
        self.MTOM_entry.insert(0, values.get('MTOM', ''))
        self.cg_position_entry.delete(0, tk.END)
        self.cg_position_entry.insert(0, values.get('cg_position', ''))

        if 'gust_velocity' in values:
            self.gust_velocity_entry.delete(0, tk.END)
            self.gust_velocity_entry.insert(0, values.get('gust_velocity', ''))

        if 'equivalent_airspeed' in values:
            self.equivalent_airspeed_entry.delete(0, tk.END)
            self.equivalent_airspeed_entry.insert(0, values.get('equivalent_airspeed', ''))


        if 'inertia_relief_weight' in values and 'inertia_relief_positions' in values:

            for widget in self.inertia_relief_frame.winfo_children():
                widget.destroy()

            inertia_weights = values.get('inertia_relief_weight', [])
            inertia_positions = values.get('inertia_relief_positions', [])

            if inertia_weights and inertia_positions:
                self.inertia_relief.set("Yes")
                self.inertia_relief_entries = []
                for idx, (position, weight) in enumerate(zip(inertia_positions, inertia_weights)):

                    ttk.Label(self.inertia_relief_frame, text=f"Inertia Relief Position {idx + 1} (m):").grid(row=idx,
                                                                                                              column=0,
                                                                                                              padx=5,
                                                                                                              pady=5)
                    position_entry = ttk.Entry(self.inertia_relief_frame)
                    position_entry.grid(row=idx, column=1, padx=5, pady=5)
                    position_entry.insert(0, position)

                    ttk.Label(self.inertia_relief_frame, text=f"Inertia Relief Weight {idx + 1} (kg):").grid(row=idx,
                                                                                                             column=2,
                                                                                                             padx=5,
                                                                                                             pady=5)
                    weight_entry = ttk.Entry(self.inertia_relief_frame)
                    weight_entry.grid(row=idx, column=3, padx=5, pady=5)
                    weight_entry.insert(0, weight)


                    self.inertia_relief_entries.append((position_entry, weight_entry))


        self.fiber_volume_fraction_skin_entry.delete(0, tk.END)
        self.fiber_volume_fraction_skin_entry.insert(0, values.get('fiber_volume_fraction_skin', ''))

        self.fiber_volume_fraction_spars_entry.delete(0, tk.END)
        self.fiber_volume_fraction_spars_entry.insert(0, values.get('fiber_volume_fraction_spars', ''))

        self.fiber_volume_fraction_skin_box_entry.delete(0, tk.END)
        self.fiber_volume_fraction_skin_box_entry.insert(0, values.get('fiber_volume_fraction_beam', ''))

        self.layer_thickness_entry.delete(0, tk.END)
        self.layer_thickness_entry.insert(0, values.get('layer_thickness', ''))

        self.show_load_case_specific_inputs()
        self.show_inertia_relief_inputs()

    def get_results(self):

        """
        This function collects the structural analysis results
        """

        results = {
            "Parameter": [
                "Load Factor", "Total Lift", "Spar 1 Height", "Spar 2 Height",
                "Skin Area", "Skin Circumference", "Material Density (Skin)", "Material Density (Spars)",
                "Spar 1 Weight", "Spar 2 Weight", "Skin Weight", "Wing Weight", "Total Weight",
                "Young's Modulus (Skin)", "Young's Modulus (Spars)", "Shear Modulus (Skin)", "Shear Modulus (Spars)",
                "Torsion Moment", "Shear Flows", "Bending Moment", "Strain", "Stresses (Top, Bottom)", "Warping", "Bending Strain Distribution",
                "Bending Stress Distribution", "Shear Force", "Shear Flow", "Bending Deflection", 'Twist'
            ],
            "Value": [
                self.load_factor, self.total_lift, self.spar1_height, self.spar2_height,
                self.skin_area, self.skin_circumference, self.material_density_skin,
                self.material_density_spars,
                self.spar1_weight, self.spar2_weight, self.skin_weight, self.wing_weight, self.total_weight,
                self.E_x_gurt, self.E_x_spar1, self.G_xy_gurt, self.G_xy_spar1,
                self.torsion_moment, self.shear_flows, self.bending_moment, self.torsion_shear_strains, self.torsion_shear_stresses, self.nu,
                self.bending_strain_distribution, self.bending_stress_distribution_mpa,
                self.Q_z, self.shear_flows_shear, self.deflection, self.twist
            ]
        }
        return results
