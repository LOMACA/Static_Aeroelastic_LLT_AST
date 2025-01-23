###################################################################

# Airfoil Geometry Module

##################################################################

import numpy as np

class AirfoilGeometry:

    """
    This class takes care of building the three-cell cross-section geometry from airfoil geometry files
    """

    def __init__(self, filepath):
        self.filepath = filepath
        self.original_coordinates = self.read_airfoil_data(filepath)
        self.coordinates = list(self.original_coordinates)
        self.chord_length = 1.0
        self.skin_thickness = 0.001
        self.skin_thickness_box = 0.002
        self.partition1 = 0.3
        self.partition2 = 0.7
        self.apply_skin_thickness()
        self.scale_airfoil(self.chord_length)
        self.upper_surface, self.lower_surface = self.separate_surfaces()
        self.max_thickness, self.max_thickness_location = self.calculate_max_thickness()

    def read_airfoil_data(self, filepath):

        """
        This function imports airfoil geometry data
        """

        with open(filepath, 'r') as file:
            lines = file.readlines()

        coordinates = []
        for line in lines[1:]:
            x, y = map(float, line.split())
            coordinates.append((x, y))

        coordinates.append(coordinates[0])

        return coordinates

    def separate_surfaces(self):

        """
        This function separates the imported airfoil into an upper and lower surface
        """

        min_x_index = np.argmin([coord[0] for coord in self.coordinates])

        upper_surface = self.coordinates[:min_x_index + 1]

        lower_surface = self.coordinates[min_x_index:][::-1]

        return upper_surface, lower_surface

    def calculate_max_thickness(self):

        """
        This function computes the maximum thickness of the airfoil
        """

        upper_surface, lower_surface = self.upper_surface, self.lower_surface

        x_upper = np.array([coord[0] for coord in upper_surface])
        y_upper = np.array([coord[1] for coord in upper_surface])

        x_lower = np.array([coord[0] for coord in lower_surface])
        y_lower = np.array([coord[1] for coord in lower_surface])

        y_lower_interp = np.interp(x_upper, x_lower, y_lower)

        thickness = y_upper - y_lower_interp

        max_thickness = np.max(thickness)
        max_thickness_location = x_upper[np.argmax(thickness)]

        return max_thickness, max_thickness_location

    def scale_airfoil(self, new_chord_length):

        """
        This function scales the airfoil chords
        """

        scale_factor = new_chord_length

        self.original_coordinates = [(x * scale_factor, y * scale_factor) for x, y in self.original_coordinates]
        self.coordinates = [(x * scale_factor, y * scale_factor) for x, y in self.coordinates]

        self.upper_surface = [(x * scale_factor, y * scale_factor) for x, y in self.upper_surface]
        self.lower_surface = [(x * scale_factor, y * scale_factor) for x, y in self.lower_surface]

        self.chord_length = new_chord_length
        self.max_thickness, self.max_thickness_location = self.calculate_max_thickness()

    def apply_skin_thickness(self):

        """
        This function simulates a skin thickness by projecting a centroid of the original airfoil outward by a scale
        factor
        """

        perimeter = self.calculate_perimeter(self.coordinates)
        scale_factor = (perimeter - 2 * self.skin_thickness) / perimeter

        centroid = np.mean(self.coordinates, axis=0)
        self.coordinates = centroid + scale_factor * (self.coordinates - centroid)

        partition1_x = self.chord_length * self.partition1
        partition2_x = self.chord_length * self.partition2
        for i, (x, y) in enumerate(self.coordinates):
            if partition1_x <= x <= partition2_x:

                scale_factor_box = (self.skin_thickness_box - self.skin_thickness) / 2
                self.coordinates[i] = (x, y + np.sign(y) * scale_factor_box)


        self.upper_surface, self.lower_surface = self.separate_surfaces()
        self.max_thickness, self.max_thickness_location = self.calculate_max_thickness()

    def reset_airfoil(self):

        """
        This function resets the airfoil to the original dimensions
        """

        self.coordinates = np.array(self.original_coordinates)
        self.upper_surface, self.lower_surface = self.separate_surfaces()
        self.max_thickness, self.max_thickness_location = self.calculate_max_thickness()

    def calculate_perimeter(self, coordinates):

        """
        This function computes the perimeter of the airfoil coordinates
        """

        perimeter = 0
        for i in range(1, len(coordinates)):
            x1, y1 = coordinates[i - 1]
            x2, y2 = coordinates[i]
            perimeter += np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        perimeter += np.sqrt(
            (coordinates[0][0] - coordinates[-1][0]) ** 2 + (coordinates[0][1] - coordinates[-1][1]) ** 2)
        return perimeter


    def update_partitions(self, partition1, partition2):

        """
        This function places the spar webs upon shifting them along the airfoil chord
        """

        self.partition1 = partition1
        self.partition2 = partition2
