"""
Parameter sets of ArchiSimple (RootArchitecture), as scenarios' "parameters".

PEA, MAIZE, PEACH
        ArchiSimple's calibrations (Pagès et al., 2014, Ecological Modelling, Table 3: Pisum sativum, Zea mays, Prunus
        persica; values in parentheses there are temporary ones). Converted to the model's units by archisimple():
        lengths in m; EL from mm mm-1 day-1 (÷ 86400); GDs from day mm-2 (× 86400 × 1e6); the tissue density from
        g cm-3 (× 1e6 g m-3); ER from day-1 (÷ 86400); LDs from day mm-2 (g cm-3)-1 (× 86400), the life duration ∝ D²
        as in the paper (life_duration_diameter_exponent 2). MNP, the number of primary roots, as one primary (seminal)
        root and nodal (adventitious) roots emitted at ER from the start (MAIZE: 5 seminal roots, nodal ones from day 5,
        an assumption). ArchiSimple has no extender on the main roots' growth durations (main_roots_growth_extender 1)
        and its durations are GDs·d². Not in the model: the tropism G, the soil's effects. Durations are in thermal time
        at the temperature factor of the soil temperature.
WHEAT   winter wheat (cv. Soissons): the architecture's parameters of RhizoDep's calibration in Wheat-BRIDGES (its
        scenario table's column GB_soil_1.0, rhizodep's WHEAT_SOISSON) in ArchiSimple's form: the life duration ∝ D²,
        LDs such that a root of D_ini lives as long as there (56.7 days); the main roots' durations extended 100 times,
        as there (GDs·d² alone would stop a 0.7 mm primary root after about 4 days: wheat's main roots are thinner than
        maize's); the nodal roots 0.95 times D_ini; laterals up to the third order, as there.
SPECIES {name: parameters}: the four.
"""
DAY = 86400.


def archisimple(Di, EL, Dmin, GDs, IPD, RDM, CVDD, SGC, LDs, RTD, n_seminal_roots=1, n_adventitious_roots=0, ER=None,
                starting_time_for_adventitious_roots_emergence=0., emergence_delay=5., **options):
    """A scenario from ArchiSimple's parameters in Pagès et al.'s (2014) units: Di, Dmin, IPD in mm; EL in mm mm-1
    day-1; GDs in day mm-2; LDs in day mm-2 (g cm-3)-1; RTD in g cm-3; ER in day-1; the delays in days. *options*:
    other parameters of the component, in its own units (they override the scenario's)."""
    scenario = dict(D_ini=Di * 1e-3, EL=EL / DAY, Dmin=Dmin * 1e-3, GDs=GDs * DAY * 1e6, IPD=IPD * 1e-3, RMD=RDM,
                    CVDD=CVDD, SGC=SGC, LDs=LDs * DAY, root_tissue_density=RTD * 1e6, new_root_tissue_density=RTD * 1e6,
                    life_duration_diameter_exponent=2., main_roots_growth_extender=1,
                    n_seminal_roots=n_seminal_roots, n_adventitious_roots=n_adventitious_roots,
                    D_sem_to_D_ini_ratio=1., D_adv_to_D_ini_ratio=1., radial_growth=SGC > 0., random=True,
                    starting_time_for_adventitious_roots_emergence=starting_time_for_adventitious_roots_emergence * DAY,
                    emergence_delay=emergence_delay * DAY)
    if ER is not None:
        scenario["ER"] = ER / DAY
    scenario.update(options)
    return scenario


PEA = archisimple(Di=1.1, EL=18, Dmin=0.19, GDs=120, IPD=7.8, RDM=0.50, CVDD=0.18, SGC=0.9, LDs=3000, RTD=0.05)
MAIZE = archisimple(Di=1.2, EL=51, Dmin=0.14, GDs=50, IPD=2.0, RDM=0.12, CVDD=0.30, SGC=0.0, LDs=3000, RTD=0.08,
                    n_seminal_roots=5, n_adventitious_roots=35, ER=0.8,
                    starting_time_for_adventitious_roots_emergence=5.)
PEACH = archisimple(Di=1.2, EL=14, Dmin=0.15, GDs=150, IPD=3.1, RDM=0.27, CVDD=0.27, SGC=0.9, LDs=3000, RTD=0.15,
                    n_seminal_roots=1, n_adventitious_roots=9, ER=0.4)

PEA = archisimple(Di=1.1, EL=18, Dmin=0.19, GDs=120, IPD=7.8, RDM=0.50, CVDD=0.18, SGC=0.9, LDs=3000, RTD=0.05)
MAIZE = archisimple(Di=1.2, EL=51, Dmin=0.14, GDs=50, IPD=2.0, RDM=0.12, CVDD=0.30, SGC=0.0, LDs=3000, RTD=0.08,
                    n_seminal_roots=5, n_adventitious_roots=35, ER=0.8,
                    starting_time_for_adventitious_roots_emergence=5.)
PEACH = archisimple(Di=1.2, EL=14, Dmin=0.15, GDs=150, IPD=3.1, RDM=0.27, CVDD=0.27, SGC=0.9, LDs=3000, RTD=0.15,
                    n_seminal_roots=1, n_adventitious_roots=9, ER=0.4)
WHEAT = archisimple(Di=0.7, EL=56.17, Dmin=0.122, GDs=8.6282, IPD=4.74, RDM=0.57, CVDD=0.16, SGC=0.0, LDs=826.72,
                    RTD=0.14, n_seminal_roots=5, n_adventitious_roots=50, ER=1.2,
                    starting_time_for_adventitious_roots_emergence=16.8274, emergence_delay=3.27,
                    main_roots_growth_extender=100, D_adv_to_D_ini_ratio=0.95, initial_segment_length=0.001,
                    initial_apex_length=0., segment_length=0.0045, root_order_limitation=True,
                    root_order_treshold=3)
SPECIES = {"wheat (Soissons)": WHEAT, "pea": PEA, "maize": MAIZE, "peach": PEACH}
