import unittest
import scripts.constants as constants
import scripts.material as material
import scripts.particle as particle
import matplotlib.pyplot as pyplot

class TestMaterial(unittest.TestCase):

    def get_sigma(self, energy, pid, mat_name):
        mat = material.Material()
        mat.set_material(mat_name, 10.)
        ion = particle.Particle.new_from_ke(energy, pid)
        sigma = mat.stripping_cross_section_nakai(ion)
        return sigma

    def test_stripping_cross_section_1(self):
        pid = constants.get_pid("H-")
        for energy, ref in (1e-3, 6e-16), (1.0, 4e-17):
            sigma = self.get_sigma(energy, pid, "gaseous_helium")
            print "E:", energy, "sigma:", sigma, "ref sigma:", ref, "fractional error", abs(1-sigma/ref), 0.2
            self.assertLess(abs(1-sigma/ref), 0.2)
        self.plot()

    def test_stripping_cross_section_2(self):
        pid = constants.get_pid("H")
        for energy, ref in [(1e-3, 4e-17),
                            (1.0, 1.5e-17),
                            (3.0, 4.7e-18),
                            (11.0, 1.2e-18)]:
            sigma = self.get_sigma(energy, pid, "gaseous_helium")
            print "E:", energy, "sigma:", sigma, "ref sigma:", ref, "fractional error", abs(1-sigma/ref), 0.2
            self.assertLess(abs(1-sigma/ref), 0.2)

    def test_strip(self):
        helium = material.Material()
        helium.set_material("gaseous_helium", 10.)
        aluminium = material.Material()
        aluminium.set_material("aluminium", 10.)
        carbon = material.Material()
        carbon.set_material("carbon", 10.)
        carbon_dioxide = material.Material()
        carbon_dioxide.set_material("carbon_dioxide", 10.)
        nitrogen = material.Material()
        nitrogen.set_material("gaseous_nitrogen", 10.)
        for energy, algorithm in (181, "saha"), (3, "nakai"):
            print "At", energy, "MeV"
            h_minus   = particle.Particle.new_from_ke(energy, constants.get_pid("H-"))
            h_neutral = particle.Particle.new_from_ke(energy, constants.get_pid("H"))
            for foil, thickness in [(helium, 1.),
                                    (aluminium, 12*1e-4),
                                    (carbon, 48*1e-6/carbon.density),
                                    (carbon_dioxide, 1.),
                                    (nitrogen, 1.e-4)]:
                foil.stripping_algorithm = algorithm
                print "For", thickness, "cm", foil.name, "using", foil.stripping_algorithm, ":"
                print "  Sigma_-1,0:              ", foil.stripping_cross_section(h_minus)
                print "  Sigma_0,1:               ", foil.stripping_cross_section(h_neutral)
                print "  Probability of H- -> H0: ", foil.strip(h_minus, thickness)
                print "  Probability of H -> H+:  ", foil.strip(h_neutral, thickness)
                print "  Energy loss:             ", foil.energy_loss_dz(h_minus)*thickness
            print
 
    def test_energy_loss(self):
        lh2 = material.Material()
        lh2.set_material("liquid_hydrogen", 10.)
        for p, ref_de in [(1.527E+02, 4.870),
                          (1.994E+02, 4.385),
                          (2.218E+02, 4.267)]:
            #ref_de is in  MeV cm^2/g
            muon = particle.Particle.new_from_momentum(p, -13)
            test_de = lh2.energy_loss_dz(muon)/lh2.density
            print lh2.name, round(muon.get_kinetic_energy()), ref_de, test_de

    def test_zzz_make_plots(self):
        if __name__ != "__main__":
            return
        print "Close plots to continue"
        pyplot.show()

    matplot_objects = []


if __name__ == "__main__":
    unittest.main()
