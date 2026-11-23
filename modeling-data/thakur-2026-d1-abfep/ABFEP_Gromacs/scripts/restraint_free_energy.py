import math

#===================================================================================================
# INPUTS
#===================================================================================================

K = 8.314472 * 0.001  # Gas constant in kJ/mol/K
V = 1.66              # Standard volume in nm^3

T      = 310.0        # Temperature in Kelvin
r0     = 0.545        # Distance in nm
thA    = 90.0         # Angle in degrees
thB    = 65.6        # Angle in degrees

K_r    = 418.4       # Force constant for distance (kJ/mol/nm^2)
K_thA  = 209.2        # Force constant for angle (kJ/mol/rad^2)
K_thB  = 209.2
K_phiA = 209.2        # Force constant for dihedral (kJ/mol/rad^2)
K_phiB = 209.2
K_phiC = 209.2

#===================================================================================================
# BORESCH FORMULA
#===================================================================================================

# Convert angles from degrees to radians
thA = math.radians(thA)
thB = math.radians(thB)

arg = (
    (8.0 * math.pi**2.0 * V) / (r0**2.0 * math.sin(thA) * math.sin(thB))
    * (
        ((K_r * K_thA * K_thB * K_phiA * K_phiB * K_phiC) ** 0.5)
        / ((2.0 * math.pi * K * T) ** 3.0)
    )
)

dG = -K * T * math.log(arg)

print("dG_off = %8.3f kcal/mol" % (dG / 4.184))
print("dG_on  = %8.3f kcal/mol" % (-dG / 4.184))


