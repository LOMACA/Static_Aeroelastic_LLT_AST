########################################################################################

# Coupled static aeroelastic model (two-way coupling)

########################################################################################

import numpy as np
import tkinter as tk
from tkinter import ttk
import logging
from llt import LiftingLineTheory

class AeroElasticModel(ttk.Frame):

    """
    This class takes care of coupling the LLT aerodynamic and analytic structural solvers
    in two-way coupling to obtain a static aeroelastic tool
    """

    def __init__(self, parent, isa_display, llt_display, airfoil_data_display, airfoil_geometry_display, wing_geometry, composite_model, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        self.isa_display = isa_display
        self.llt_display = llt_display
        self.airfoil_data_display = airfoil_data_display
        self.airfoil_geometry_display = airfoil_geometry_display
        self.wing_geometry = wing_geometry
        self.composite_model = composite_model

        self.create_input_widgets()
        self.create_result_display()

    def create_input_widgets(self):

        """
        This function creates input fields for the static aeroelastic analysis
        """

        input_frame = ttk.LabelFrame(self, text="Inputs")
        input_frame.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")

        ttk.Label(input_frame, text="Initial Angle of Attack (degrees):").grid(row=0, column=0, padx=5, pady=5)
        self.alpha_entry = ttk.Entry(input_frame)
        self.alpha_entry.grid(row=0, column=1, padx=5, pady=5)

        ttk.Button(input_frame, text="Compute Aeroelastic Coupling", command=lambda: self.compute_aeroelastic_coupling(float(self.alpha_entry.get()))).grid(row=1, column=0, columnspan=2, padx=5, pady=5)

    def create_result_display(self):

        """
        This function creates a result window for the convergence information
        """

        self.result_frame = ttk.LabelFrame(self, text="Results")
        self.result_frame.grid(row=0, column=1, padx=5, pady=5, sticky="nsew")

        self.result_scrollbar = ttk.Scrollbar(self.result_frame)
        self.result_scrollbar.grid(row=0, column=1, sticky="ns")

        self.result_text = tk.Text(self.result_frame, wrap=tk.WORD, width=60, height=20)
        self.result_text.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")

        self.result_text.config(yscrollcommand=self.result_scrollbar.set)
        self.result_scrollbar.config(command=self.result_text.yview)

    def normalize_deflections(self, deflections):

        """
        This function normalizes the bending deflection by the chord
        """

        span_positions = np.arange(len(deflections))
        chord_lengths = np.array([self.wing_geometry.get_chord_length(pos) for pos in span_positions])
        normalized_deflections = deflections / chord_lengths
        return normalized_deflections

    def compute_aeroelastic_coupling(self, initial_alpha, max_iterations=1000, tolerance=2000):

        """
        This function solves the coupled aeroelastic model
        """

        self.result_text.delete("1.0", tk.END)
        results_text = ""

        logging.basicConfig(level=logging.INFO, format='%(message)s')
        logger = logging.getLogger()

        current_alpha = initial_alpha

        relaxation_factor = 0.1
        deflection_adjustment_weight = 0.2
        twist_adjustment_weight = 0.2

        mtom = self.composite_model.get_mtom() * 9.81
        load_factor = self.composite_model.load_factor
        target_lift = load_factor * mtom

        max_alpha = 20
        min_alpha = -20
        max_alpha_change = 3

        for iteration in range(max_iterations):

            logger.info(f"\nIteration {iteration + 1} - LLT Inputs:")

            lifting_line = LiftingLineTheory(
                v_inf=self.llt_display.v_inf,
                c_root=self.llt_display.c_root,
                c_tip=self.llt_display.c_tip,
                b=self.llt_display.b,
                n=self.llt_display.n,
                cl0=self.airfoil_data_display.root_cl0,
                cd0=self.airfoil_data_display.root_cd0,
                a0=self.llt_display.a0,
                alpha0=self.llt_display.alpha0,
                rho=self.isa_display.rho
            )

            llt_results = lifting_line.compute([current_alpha])
            lift_distribution = llt_results[0]['L_dist']
            L_total = np.sum(llt_results['Lift'])

            self.composite_model.compute(lift_distribution)
            deflections = self.composite_model.deflection
            twists = self.composite_model.twist

            normalized_deflections = self.normalize_deflections(deflections)

            max_deflection = 4
            max_twist = 5
            normalized_deflections = np.clip(normalized_deflections, -max_deflection, max_deflection)
            twists = np.clip(twists, -max_twist, max_twist)

            lift_weight_diff = L_total - target_lift
            logger.info(f"Iteration {iteration + 1}: Lift-Weight Difference = {lift_weight_diff:.4f} N")

            if abs(lift_weight_diff) < tolerance:
                logger.info(f"Converged after {iteration + 1} iterations.")
                break

            d_alpha_deflection_adjustment = deflection_adjustment_weight * np.degrees(np.mean(np.gradient(normalized_deflections)))
            d_alpha_twist_adjustment = twist_adjustment_weight * np.degrees(np.mean(twists))

            alpha_change = relaxation_factor * (d_alpha_deflection_adjustment + d_alpha_twist_adjustment)
            alpha_change = np.clip(alpha_change, -max_alpha_change, max_alpha_change)

            current_alpha += alpha_change
            current_alpha = np.clip(current_alpha, min_alpha, max_alpha)

            results_text += f"Iteration {iteration + 1}:\n"
            results_text += f"Angle of Attack: {current_alpha:.4f} degrees\n"
            results_text += f"Total Lift: {L_total:.2f} N\n"
            results_text += f"Target Lift (n*W): {target_lift:.2f} N\n"
            results_text += f"Twist: {twists} degrees\n"
            results_text += f"Deflection: {deflections} m\n\n"

        else:
            results_text += f"Failed to converge within {max_iterations} iterations.\n"

        self.result_text.insert(tk.END, results_text)

    def get_input_values(self):
        pass

    def set_input_values(self, values):
        pass

    def get_results(self):
        pass
