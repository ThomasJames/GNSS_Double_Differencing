from math import sqrt


class PartialDiffCalc:
    def __init__(self, ref_station, corresponding_sat, sat_ref, wavelength):
        self.X_3A = ref_station[0]
        self.Y_3A = ref_station[1]
        self.Z_3A = ref_station[2]

        self.X_s = corresponding_sat[0]
        self.Y_s = corresponding_sat[1]
        self.Z_s = corresponding_sat[2]

        self.X_s_ref = sat_ref[0]
        self.Y_s_ref = sat_ref[1]
        self.Z_s_ref = sat_ref[2]

        self.wavelength = wavelength

    def x_diff(self):
        result = 1 / self.wavelength * \
                 (
                         (self.X_3A - self.X_s) /
                         (sqrt((self.X_s - self.X_3A) ** 2 +
                               (self.Y_s - self.Y_3A) ** 2 +
                               (self.Z_s - self.Z_3A) ** 2))
                         -
                         (self.X_3A - self.X_s_ref) /
                         (sqrt(self.X_s_ref - self.X_3A) ** 2 +
                          (self.Y_s_ref - self.Y_3A) ** 2 +
                          (self.Z_s_ref - self.Z_3A) ** 2)
                 )
        return float(result)

    def y_diff(self):
        result = 1 / self.wavelength * \
                 (
                         (self.Y_3A - self.Y_s) /
                         (sqrt((self.X_s - self.X_3A) ** 2 +
                               (self.Y_s - self.Y_3A) ** 2 +
                               (self.Z_s - self.Z_3A) ** 2))
                         -
                         (self.Y_3A - self.Y_s_ref) /
                         (sqrt(self.X_s_ref - self.X_3A) ** 2 +
                          (self.Y_s_ref - self.Y_3A) ** 2 +
                          (self.Z_s_ref - self.Z_3A) ** 2)
                 )
        return float(result)

    def z_diff(self):
        result = 1 / self.wavelength * \
                 (
                         (self.Z_3A - self.Z_s) /
                         (sqrt((self.X_s - self.X_3A) ** 2 +
                               (self.Y_s - self.Y_3A) ** 2 +
                               (self.Z_s - self.Z_3A) ** 2))
                         -
                         (self.Z_3A - self.Z_s_ref) /
                         (sqrt(self.X_s_ref - self.X_3A) ** 2 +
                          (self.Y_s_ref - self.Y_3A) ** 2 +
                          (self.Z_s_ref - self.Z_3A) ** 2)
                 )
        return float(result)
