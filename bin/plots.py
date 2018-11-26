import scripts.constants as constants
import scripts.material as material
import scripts.particle as particle
import matplotlib.pyplot as pyplot

########### CROSS SECTIONS ################
def get_sigma(energy, pid, mat_name, algorithm):
    mat = material.Material()
    mat.set_material(mat_name, 10.)
    mat.stripping_algorithm = algorithm
    ion = particle.Particle.new_from_ke(energy, pid)
    sigma = mat.stripping_cross_section(ion)
    return sigma

def plot_sigma(ion, plot_id, fig):
    pid = constants.get_pid(ion)
    ax1 = fig.add_subplot(plot_id)
    for mat_name, algo in [("gaseous_helium", "nakai"),
                           ("gaseous_nitrogen", "nakai"),
                           ("carbon_dioxide", "nakai"),
                           ("carbon", "saha")]:
        labels = {"nakai":"N-S", "saha":"Saha et al."}
        energy_list = sorted([10**(0.25*i)*1e-6 for i in range(13, 37)]+[3.])
        sigma_list = [get_sigma(energy, pid, mat_name, algo) for energy in energy_list]
        a_plot = ax1.plot(energy_list, sigma_list, label=mat_name+" - "+labels[algo])
        ref_sigma = get_sigma(3., pid, mat_name, algo)
        if ion == "H-":
            sub = "_{-10}"
        else:
            sub = "_{01}"
        ax1.plot([min(energy_list), max(energy_list)],
                 [ref_sigma, ref_sigma],
                 label="$\sigma"+sub+"(3 MeV)=$"+format(ref_sigma, '.2e'),
                 linewidth=1,
                 linestyle='dashed',
                 color=a_plot[0].get_color())
        ax1.plot([3., 3.],
                 [min(sigma_list), max(sigma_list)],
                 label=None,
                 linewidth=1,
                 linestyle='dashed',
                 color=a_plot[0].get_color())
    ax1.set_xscale('log')
    ax1.set_yscale('log')

def labels(fig):
    ax1 = fig.add_subplot(221)
    ax1.text(0.5, 0.9, 'Stripping of H$^{-}$', transform = ax1.transAxes)
    ax1.set_ylabel('Cross section [cm$^{2}$]')
    ax1.yaxis.set_ticks([10.**i for i in range(-20, -11)])
    ax1.legend(bbox_to_anchor=(1.01, 1.1, 1.49, 0.2))
    ax2 = fig.add_subplot(223)
    ax2.text(0.5, 0.9, 'Stripping of H$^{0}$', transform = ax2.transAxes)
    ax2.yaxis.set_ticks([10.**i for i in range(-20, -11)])
    ax2.set_ylabel('Cross section [cm$^{2}$]')
    ax2.set_xlabel('Kinetic Energy [MeV]')
    ax2.legend(bbox_to_anchor=(1.01, 1.0, 1.49, 0.1))

def cross_sections():
    fig = pyplot.figure(figsize=(6, 4), dpi=250)
    plot_sigma("H-", 221, fig)
    plot_sigma("H", 223, fig)
    labels(fig)
    fig.savefig("plots/cross_section.png")

########### ENERGY LOSS ##############


def get_dedz(mat_name, energy):
    mat = material.Material()
    mat.set_material(mat_name, 10.)
    ion = particle.Particle.new_from_ke(energy, constants.get_pid("H-"))
    dedz = mat.energy_loss_dz(ion)
    return dedz

def energy_loss():
    fig = pyplot.figure()
    ax1 = fig.add_subplot(111)
    for mat_name in [("gaseous_helium"),
                     ("gaseous_nitrogen"),
                     ("carbon_dioxide"),
                     ("carbon")]:
        energy_list = [10**(0.25*i)*1e-6 for i in range(24, 37)]
        dedz_list = [-get_dedz(mat_name, energy) for energy in energy_list]
        ax1.plot(energy_list, dedz_list, label=mat_name)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.legend() #bbox_to_anchor=(0., 0.1, 0.5, 0.5))
    ax1.set_xlabel('Kinetic Energy [MeV]')
    ax1.set_ylabel('dE/dx [MeV/cm]')
    fig.savefig("plots/dedx.png")

########### STRIPPING IN BEAM PIPE ###########
def lifetime(distance_factor):
    # distance in cm
    fig = pyplot.figure()
    ax1 = fig.add_subplot(111)
    mat = material.Material()
    mat.set_material("gaseous_nitrogen", 10.)
    mat.stripping_algorithm = "nakai"
    energy_list = [3., 10., 20., 30.]
    for energy in energy_list:
        ion = particle.Particle.new_from_ke(energy, constants.get_pid("H-"))
        pressure_list = [10**(i/16.-6) for i in range(80)]
        survival_list = [mat.strip(ion, pressure) for pressure in pressure_list]
        pressure_list = [pressure/distance_factor for pressure in pressure_list]
        ax1.plot(pressure_list, survival_list, label=str(energy)+" MeV")
    ax1.set_xscale('log')
    #ax1.set_yscale('log')
    ax1.legend() #bbox_to_anchor=(0., 0.1, 0.5, 0.5))
    ax1.set_title('H- ion survival in gaseous nitrogen')
    ax1.set_xlabel('Pressure over '+str(distance_factor*1e-2)+' m flight path [bar]')
    ax1.set_ylabel('Ionisation probability')
    fig.savefig("plots/lifetime.png")

########## FOIL STRIPPING THICKNESS ##########
def get_strip(target_fraction, tolerance, function):
    d0, f0 = 1e-16, 0.
    d1, f1 = 1., 1.
    while f1-f0 > tolerance:
        dtest = (d1+d0)/2.
        ftest = function(dtest)
        if ftest < target_fraction:
            d0 = dtest
            f0 = ftest
        else:
            d1 = dtest
            f1 = ftest
        #print d0, d1, "**", f0, f1
    return (d1+d0)/2., (f1+f0)/2.

def stripping():
    print 'Species'.ljust(8), "Material".ljust(20), '"Thickness to strip 99.9% [mm]"'
    for species in ["H-", "H"]:
        for mat_name, algo in [("gaseous_helium", "nakai"),
                              ("gaseous_nitrogen", "nakai"),
                              ("carbon", "saha")]:
            mat = material.Material()
            mat.set_material(mat_name, 10.)
            mat.stripping_algorithm = algo
            ion = particle.Particle.new_from_ke(3, constants.get_pid(species))
            thickness, stripping = get_strip(0.999, 0.00001, lambda dx: mat.strip(ion, dx))
            thickness *= 10. # convert cm -> mm
            print species.ljust(8), mat_name.ljust(20), format(thickness, '.2e')


def foil_de():
    print "Material".ljust(20), '"Thickness to lose 10% [mm]"'
    for mat_name, algo in [("gaseous_helium", "nakai"),
                          ("gaseous_nitrogen", "nakai"),
                          ("carbon", "saha"),
                          ("aluminium", "saha")]:
        mat = material.Material()
        mat.set_material(mat_name, 10.)
        mat.stripping_algorithm = algo
        ion = particle.Particle.new_from_ke(3, constants.get_pid("H-"))
        dedz = mat.energy_loss_dz(ion)
        thickness = -0.3/dedz*10.
        print mat_name.ljust(20), format(thickness, '.2e')


def main():
    cross_sections()
    energy_loss()
    lifetime(1.)
    lifetime(30.*100.*25000)
    stripping()
    foil_de()
    pyplot.show()

if __name__ == "__main__":
    main()