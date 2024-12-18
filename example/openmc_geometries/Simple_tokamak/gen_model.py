import openmc
import openmc.model
import os
import numpy as np
from openmc.stats import (
    Discrete,
    Uniform,
    CylindricalIndependent,
    Isotropic,
    muir
)

def SimpleTokamak(radius = 500, first_wall_thicknesses = 5, blanket_thicknesses=100, center_column_thicknesses=50, shield_thickness = 100):
    """
    This function creates an sphere with two layers.
    Returns Tritium Breeding Ratio and Neutron Leakage.
    The inside sphere is void with 1cm of radius.
    Next layer is Tungsten.
    Last layer is PbLi.
    Variable Radius is thickness of the Tungsten Layer
    """
#     radius = radius + 1
#     Material 0 is W
#     Material 1 is LiPB
#     Material 2 is Fe
#     Material 3 is Copper


    #MATERIALS#

    mats = openmc.Materials()

    tungsten = openmc.Material(name='Tungsten', material_id=1)
    tungsten.set_density('g/cm3', 19.0)
    tungsten.add_element('W', 1.0)
    mats.append(tungsten)

    enrichment_fraction = 0.90
    breeder_material = openmc.Material(name='PbLi', material_id=2) #Pb84.2Li15.8 with enrichment of Li6
    breeder_material.set_density('g/cm3', 9.956685348)
    breeder_material.add_element('Pb', 84.2,'ao')
    breeder_material.add_nuclide('Li6', enrichment_fraction*15.8, 'ao')
    breeder_material.add_nuclide('Li7', (1.0-enrichment_fraction)*15.8, 'ao')
    mats.append(breeder_material)

    eurofer = openmc.Material(name='EUROFER97',material_id=3)
    eurofer.set_density('g/cm3', 7.75)
    eurofer.add_element('Fe', 89.067, percent_type='wo')
    eurofer.add_element('C', 0.11, percent_type='wo')
    eurofer.add_element('Mn', 0.4, percent_type='wo')
    eurofer.add_element('Cr', 9.0, percent_type='wo')
    eurofer.add_element('Ta', 0.12, percent_type='wo')
    eurofer.add_element('W', 1.1, percent_type='wo')
    eurofer.add_element('N', 0.003, percent_type='wo')
    eurofer.add_element('V', 0.2, percent_type='wo')
    mats.append(eurofer)

    copper = openmc.Material(name='coppper',material_id=4)
    copper.set_density('g/cm3', 8.96)
    copper.add_element('Cu', 1.0, percent_type='ao')

    mats.append(copper)
    mats.export_to_xml()

    #GEOMETRY#

    surface1 = openmc.Sphere(r=radius)
    surface2 = openmc.Sphere(r=radius + first_wall_thicknesses)
    surface3 = openmc.Sphere(r=radius + first_wall_thicknesses + blanket_thicknesses)
    surface4 = openmc.Sphere(r=radius + first_wall_thicknesses + blanket_thicknesses + shield_thickness, boundary_type="vacuum")
    surface5 = openmc.ZCylinder(r=center_column_thicknesses)
    
    # region1 = -surface1 & +surface5              # plasma
    # region2 = +surface1 & -surface2 & +surface5  # first wall
    # region3 = +surface2 & -surface3 & +surface5  # blanket
    # region4 = +surface3 & -surface4 & +surface5  # shield
    # region5 = -surface4 & -surface5              # center column

    region1 = -surface1 & +surface5              # plasma
    region2 = +surface1 & -surface2              # first wall
    region3 = +surface2 & -surface3              # blanket
    region4 = +surface3 & -surface4              # shield
    region5 = -surface4 & -surface5              # center column

    cell1 = openmc.Cell(region=region1)
    cell2 = openmc.Cell(region=region2)
    cell2.fill = tungsten
    cell3 = openmc.Cell(region=region3)
    cell3.fill = breeder_material
    cell4 = openmc.Cell(region=region4)
    cell4.fill = eurofer
    cell5 = openmc.Cell(region=region5)
    cell5.fill = copper

    geometry = openmc.Geometry([cell1, cell2, cell3, cell4, cell5])
    geometry.export_to_xml()


    #SETTINGS#

    batches = 10
    inactive = 0
    particles = 100_000

    my_source = openmc.Source()

    ring_rad = center_column_thicknesses + (radius-center_column_thicknesses)/2
    source_radius = Discrete([ring_rad], [1])
    z_values = Discrete([0], [1])

    angle = Uniform(a=0.0, b= 2 * 3.14159265359)
    my_source.space = CylindricalIndependent(
        r=source_radius, phi=angle, z=z_values, origin=(0.0, 0.0, 0.0)
    )
    my_source.angle = Isotropic()
    my_source.energy = muir(e0=14080000.0, m_rat=5.0, kt=20000.0)

    sett = openmc.Settings()
    sett.batches = batches
    sett.inactive = inactive
    sett.particles = particles
    # sett.output = {'tallies': False}
    sett.run_mode = 'fixed source'
    sett.source = my_source
    sett.export_to_xml()

    #TALLIES#

    tallies = openmc.Tallies()

    cell_filter = openmc.CellFilter(cell3)
    tbr_tally = openmc.Tally(name='TBR')
    tbr_tally.filters = [cell_filter]
    tbr_tally.scores = ['(n,Xt)'] # MT 205 is the (n,Xt) reaction where X is a wildcard, if MT 105 or (n,t) then some tritium production will be missed, for example (n,nt) which happens in Li7 would be missed
    tallies.append(tbr_tally)

    filter = openmc.SurfaceFilter(surface4)
    leakage_tally = openmc.Tally(name='leakage')
    leakage_tally.filters = [filter]
    leakage_tally.scores = ['current']
    tallies.append(leakage_tally)

    total_rad = radius + first_wall_thicknesses + blanket_thicknesses + shield_thickness

    tally_flux = openmc.Tally(name='Flux_mesh')
    mesh = openmc.RegularMesh()
    mesh.lower_left = (-total_rad, -total_rad, -total_rad)
    mesh.upper_right = (total_rad, total_rad, total_rad)

    mesh.dimension = (200, 200, 200)
    # tally_flux.estimator = 'collision'
    # tally_heating.estimator = 'collision'
    tally_flux.scores = ['flux']

    particle_filter = openmc.ParticleFilter(['neutron'])

    mesh_filter = openmc.MeshFilter(mesh)
    tally_flux.filters = [particle_filter, mesh_filter]
    tallies.append(tally_flux)

    tbr_tally_mesh = openmc.Tally(name='TBR_mesh')
    tbr_tally_mesh.filters = [cell_filter, mesh_filter]
    tbr_tally_mesh.scores = ['(n,Xt)'] # MT 205 is the (n,Xt) reaction where X is a wildcard, if MT 105 or (n,t) then some tritium production will be missed, for example (n,nt) which happens in Li7 would be missed
    tallies.append(tbr_tally_mesh)

    tallies.export_to_xml()

    model = openmc.model.Model(geometry, mats, sett, tallies)

    #RUN#
    # model.run(output=False)
    return model


SimpleTokamak()