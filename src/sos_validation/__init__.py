"""Validate adjeff's six 5S radiative terms against the SOS code.

Modules
-------
    optics
        The CAMS aerosol optics both codes read.
    atmosphere
        The optical depths of the Smart-G atmosphere adjeff builds.
    sos
        Run SOS and parse its outputs.
    five_s
        The six terms from SOS runs over several albedos.
    adjeff_side
        The six terms from adjeff's ``RadiativePipeline``.
    validate
        The sweeps, the comparison and the command line.
    report
        The Markdown report.
    plots
        The figures.
"""
