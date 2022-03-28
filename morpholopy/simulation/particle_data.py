import numpy as np
from .unitilies.helper_functions import cosmic_time_approx_Gyr
import unyt
import swiftsimio as sw
import velociraptor


class GasParticleData:
    """
    General class containing gas particle properties
    """

    def __init__(self, sim_info, halo_id):

        halo_index = sim_info.halo_data.halo_ids == halo_id

        mask = sw.mask(f"{sim_info.directory}/{sim_info.snapshot_name}")

        # The full metadata object is available from within the mask
        size = unyt.unyt_array([0.5, 0.5, 0.5], 'Mpc')

        x = sim_info.halo_data.xminpot[halo_index]
        y = sim_info.halo_data.yminpot[halo_index]
        z = sim_info.halo_data.zminpot[halo_index]

        origin = unyt.unyt_array([x, y, z], 'Mpc') / sim_info.a  # to comoving

        # region is a 3x2 list [[left, right], [bottom, top], [front, back]]
        region = [[-0.5 * b + o, 0.5 * b + o] for b, o in zip(size, origin)]

        # Constrain the mask
        mask.constrain_spatial(region)

        # Now load the snapshot with this mask
        data = sw.load(f"{sim_info.directory}/{sim_info.snapshot_name}", mask=mask)

        gas_ids = data.stars.particle_ids.value
        mask_gas = sim_info.select_bound_particles(halo_id, gas_ids)

        # mask_gas = sim_info.make_mask_gas(halo_id=halo_id)

        self.mass = data.gas.masses[mask_gas].value * sim_info.to_Msun_units
        #self.mass = sim_info.snapshot.gas.masses[mask_gas].value * sim_info.to_Msun_units
        self.n_parts = len(self.mass)

        # self.xcoordinate = sim_info.snapshot.gas.coordinates[mask_gas, 0].value * sim_info.a * sim_info.to_kpc_units
        # self.ycoordinate = sim_info.snapshot.gas.coordinates[mask_gas, 1].value * sim_info.a * sim_info.to_kpc_units
        # self.zcoordinate = sim_info.snapshot.gas.coordinates[mask_gas, 2].value * sim_info.a * sim_info.to_kpc_units
        #
        # self.xvelocity = sim_info.snapshot.gas.velocities[mask_gas, 0].value  # km/s
        # self.yvelocity = sim_info.snapshot.gas.velocities[mask_gas, 1].value  # km/s
        # self.zvelocity = sim_info.snapshot.gas.velocities[mask_gas, 2].value  # km/s

        self.xcoordinate = data.gas.coordinates[mask_gas, 0].value * sim_info.a * sim_info.to_kpc_units
        self.ycoordinate = data.gas.coordinates[mask_gas, 1].value * sim_info.a * sim_info.to_kpc_units
        self.zcoordinate = data.gas.coordinates[mask_gas, 2].value * sim_info.a * sim_info.to_kpc_units

        self.xvelocity = data.gas.velocities[mask_gas, 0].value  # km/s
        self.yvelocity = data.gas.velocities[mask_gas, 1].value  # km/s
        self.zvelocity = data.gas.velocities[mask_gas, 2].value  # km/s

        self.smoothing_length = data.gas.smoothing_lengths[mask_gas].value * sim_info.a * sim_info.to_kpc_units
        #self.smoothing_length = sim_info.snapshot.gas.smoothing_lengths[mask_gas].value * sim_info.a * sim_info.to_kpc_units

        # XH = sim_info.snapshot.gas.element_mass_fractions.hydrogen[mask_gas].value
        # gas_HI = sim_info.snapshot.gas.species_fractions.HI[mask_gas].value
        # gas_H2 = sim_info.snapshot.gas.species_fractions.H2[mask_gas].value * 2.0

        XH = data.gas.element_mass_fractions.hydrogen[mask_gas].value
        gas_HI = data.gas.species_fractions.HI[mask_gas].value
        gas_H2 = data.gas.species_fractions.H2[mask_gas].value * 2.0

        self.HI_mass = gas_HI * XH * self.mass
        self.H2_mass = gas_H2 * XH * self.mass

        # self.star_formation_rates = sim_info.snapshot.gas.star_formation_rates[mask_gas].value * sim_info.to_Msun_units / sim_info.to_yr_units
        # self.densities = sim_info.snapshot.gas.densities[mask_gas].value * (sim_info.a * sim_info.to_Msun_units / sim_info.to_kpc_units) ** 3
        # self.metal_mass_fractions = sim_info.snapshot.gas.metal_mass_fractions[mask_gas].value / sim_info.Zsolar

        self.star_formation_rates = data.gas.star_formation_rates[mask_gas].value * sim_info.to_Msun_units / sim_info.to_yr_units
        self.densities = data.gas.densities[mask_gas].value * (sim_info.a * sim_info.to_Msun_units / sim_info.to_kpc_units) ** 3
        self.metal_mass_fractions = data.gas.metal_mass_fractions[mask_gas].value / sim_info.Zsolar


class StarParticleData:
    """
    General class containing gas particle properties
    """

    def __init__(self, sim_info, halo_id):

        halo_index = sim_info.halo_data.halo_ids == halo_id

        mask = sw.mask(f"{sim_info.directory}/{sim_info.snapshot_name}")

        # The full metadata object is available from within the mask
        size = unyt.unyt_array([0.5, 0.5, 0.5], 'Mpc')

        x = sim_info.halo_data.xminpot[halo_index]
        y = sim_info.halo_data.yminpot[halo_index]
        z = sim_info.halo_data.zminpot[halo_index]

        origin = unyt.unyt_array([x, y, z], 'Mpc') / sim_info.a  # to comoving

        # region is a 3x2 list [[left, right], [bottom, top], [front, back]]
        region = [[-0.5 * b + o, 0.5 * b + o] for b, o in zip(size, origin)]

        # Constrain the mask
        mask.constrain_spatial(region)

        # Now load the snapshot with this mask
        data = sw.load(f"{sim_info.directory}/{sim_info.snapshot_name}", mask=mask)

        stars_ids = data.stars.particle_ids.value
        mask_stars = sim_info.select_bound_particles(halo_id, stars_ids)

        #mask_stars = sim_info.make_mask_stars(halo_id=halo_id)

        # self.mass = sim_info.snapshot.stars.masses[mask_stars].value * sim_info.to_Msun_units
        # stars_birthz = (
        #         1.0 / sim_info.snapshot.stars.birth_scale_factors[mask_stars].value - 1.0
        # )

        self.mass = data.stars.masses[mask_stars].value * sim_info.to_Msun_units
        stars_birthz = (
                1.0 / data.stars.birth_scale_factors[mask_stars].value - 1.0
        )

        if len(stars_birthz) > 1:
            self.age = cosmic_time_approx_Gyr(
                z=0.0, Omega_L=sim_info.Omega_l, Hubble_time=sim_info.hubble_time_Gyr
            ) - cosmic_time_approx_Gyr(
                z=stars_birthz, Omega_L=sim_info.Omega_l, Hubble_time=sim_info.hubble_time_Gyr
            )
        else:
            self.age = 0.0

        self.metal_mass_fractions = data.stars.metal_mass_fractions[mask_stars].value
        self.initmass = data.stars.initial_masses[mask_stars].value * sim_info.to_Msun_units
        # self.metal_mass_fractions = sim_info.snapshot.stars.metal_mass_fractions[mask_stars].value
        # self.initmass = sim_info.snapshot.stars.initial_masses[mask_stars].value * sim_info.to_Msun_units

        self.n_parts = len(self.mass)

        # self.xcoordinate = sim_info.snapshot.stars.coordinates[mask_stars, 0].value * sim_info.a * sim_info.to_kpc_units
        # self.ycoordinate = sim_info.snapshot.stars.coordinates[mask_stars, 1].value * sim_info.a * sim_info.to_kpc_units
        # self.zcoordinate = sim_info.snapshot.stars.coordinates[mask_stars, 2].value * sim_info.a * sim_info.to_kpc_units
        #
        # self.xvelocity = sim_info.snapshot.stars.velocities[mask_stars, 0].value  # km/s
        # self.yvelocity = sim_info.snapshot.stars.velocities[mask_stars, 1].value  # km/s
        # self.zvelocity = sim_info.snapshot.stars.velocities[mask_stars, 2].value  # km/s

        self.xcoordinate = data.stars.coordinates[mask_stars, 0].value * sim_info.a * sim_info.to_kpc_units
        self.ycoordinate = data.stars.coordinates[mask_stars, 1].value * sim_info.a * sim_info.to_kpc_units
        self.zcoordinate = data.stars.coordinates[mask_stars, 2].value * sim_info.a * sim_info.to_kpc_units

        self.xvelocity = data.stars.velocities[mask_stars, 0].value  # km/s
        self.yvelocity = data.stars.velocities[mask_stars, 1].value  # km/s
        self.zvelocity = data.stars.velocities[mask_stars, 2].value  # km/s

        self.baryon_max_soft = 0.5 * sim_info.baryon_max_soft * np.ones(self.n_parts)
        self.weighted_mass = self.mass * (1.2348 / self.baryon_max_soft) ** 3

        # self.oxygen = sim_info.snapshot.stars.element_mass_fractions.oxygen[mask_stars].value
        # self.iron = sim_info.snapshot.stars.element_mass_fractions.iron[mask_stars].value
        # self.magnesium = sim_info.snapshot.stars.element_mass_fractions.magnesium[mask_stars].value
        # self.hydrogen = sim_info.snapshot.stars.element_mass_fractions.hydrogen[mask_stars].value
        # self.carbon = sim_info.snapshot.stars.element_mass_fractions.carbon[mask_stars].value
        # self.silicon = sim_info.snapshot.stars.element_mass_fractions.silicon[mask_stars].value
        # self.europium = sim_info.snapshot.stars.element_mass_fractions.europium[mask_stars].value
        # self.nitrogen = sim_info.snapshot.stars.element_mass_fractions.nitrogen[mask_stars].value
        # self.neon = sim_info.snapshot.stars.element_mass_fractions.neon[mask_stars].value
        # self.iron_SNIa_fraction = sim_info.snapshot.stars.iron_mass_fractions_from_snia[mask_stars].value


        self.oxygen = data.stars.element_mass_fractions.oxygen[mask_stars].value
        self.iron = data.stars.element_mass_fractions.iron[mask_stars].value
        self.magnesium = data.stars.element_mass_fractions.magnesium[mask_stars].value
        self.hydrogen = data.stars.element_mass_fractions.hydrogen[mask_stars].value
        self.carbon = data.stars.element_mass_fractions.carbon[mask_stars].value
        self.silicon = data.stars.element_mass_fractions.silicon[mask_stars].value
        self.europium = data.stars.element_mass_fractions.europium[mask_stars].value
        self.nitrogen = data.stars.element_mass_fractions.nitrogen[mask_stars].value
        self.neon = data.stars.element_mass_fractions.neon[mask_stars].value
        self.iron_SNIa_fraction = data.stars.iron_mass_fractions_from_snia[mask_stars].value

        if (hasattr(sim_info.snapshot.stars.element_mass_fractions, 'barium')):
            # self.barium = sim_info.snapshot.stars.element_mass_fractions.barium[mask_stars].value
            # self.strontium = sim_info.snapshot.stars.element_mass_fractions.strontium[mask_stars].value
            self.barium = data.stars.element_mass_fractions.barium[mask_stars].value
            self.strontium = data.stars.element_mass_fractions.strontium[mask_stars].value
        else:
            self.barium = np.zeros(self.n_parts)
            self.strontium = np.zeros(self.n_parts)

        self.in_halo = np.zeros(self.n_parts)

        x = self.xcoordinate - sim_info.halo_data.xminpot[halo_index]
        y = self.ycoordinate - sim_info.halo_data.yminpot[halo_index]
        z = self.zcoordinate - sim_info.halo_data.zminpot[halo_index]
        r = x ** 2 + y ** 2 + z ** 2
        halo_stars = np.where(r > 10 ** 2)[0]  # further than 8kpc?
        self.in_halo[halo_stars] = np.ones(len(halo_stars))

        #self.luminosity = sim_info.snapshot.stars.luminosity.rband[mask_stars].value
