using Test, DifferentialAlgebra, Aqua

# Inspect ambiguities owned by this package; dependencies own their additional
# dispatch. Extension contracts and both fresh-process loading orders are tested
# separately. All other Aqua checks use their defaults.
Aqua.test_all(DifferentialAlgebra; ambiguities = (recursive = false,))
