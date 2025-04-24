
using Pkg

Pkg.status()
Pkg.add(Pkg.PackageSpec(; name="DrWatson", version="2.18.0"))

@quickactivate "ImmuneBoostingODEs"
Pkg.instantiate()

Pkg.status()
