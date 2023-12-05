import pprint

from utils.orbit_loader import get_all_orbits, get_orbit

orbit = get_orbit("Figure-8")
pprint.pprint(orbit)
orbit_two = get_orbit("Butterfly-I")
pprint.pprint(orbit_two)
all_orbits = get_all_orbits()
pprint.pprint(all_orbits)
