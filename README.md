# MomentumConservedExactDiagonalization.jl

A Julia package for momentum-conserved exact diagonalization of quantum many-body systems, targeting researchers in condensed matter physics.

## Features

- **Momentum-conserved exact diagonalization** with automatic block diagonalization by momentum sectors
- **Advanced analysis tools** including entanglement entropy calculations, twisted boundary conditions, and many-body Chern numbers
- **High-performance implementation** leveraging Julia's sparse matrix capabilities and parallel computation
- **Comprehensive examples** for fractional quantum Hall effect, Chern insulators, and other quantum many-body systems

## Installation

Once the package is registered in the Julia General registry, you can install it using:

```julia
using Pkg
Pkg.add("MomentumConservedExactDiagonalization")
```

For development, you can install from the local directory:

```julia
using Pkg
Pkg.add(path="path/to/MomentumConservedExactDiagonalization.jl")
```

## Quick Start

```julia
using MomentumConservedExactDiagonalization

# Define system parameters
params = EDPara(
    # Parameters will be defined as the package develops
)

# Build Hilbert space
hilbert_space = HilbertSpace(params)

# Construct Hamiltonian
hamiltonian = HamiltonianBuilder(hilbert_space, params)

# Diagonalize
results = Diagonalizer(hamiltonian)
```

## Documentation

Full documentation is available at [GitHub Pages link to be added].

## Development Status

This package is currently in development. Core functionality will be implemented in the following phases:

1. **Package Infrastructure** ✅ - Complete
2. **Core ED Engine** - In Progress
3. **Advanced Analysis Features** - Planned
4. **Example Systems** - Planned
5. **Documentation & Release** - Planned

## Contributing

Contributions are welcome! Please see the development roadmap in the project management files.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this package in your research, please cite:

```bibtex
@software{momentum_conserved_ed,
  title={MomentumConservedExactDiagonalization.jl},
  author={Zou Bo},
  year={2025},
  url={https://github.com/Zou-Bo/MomentumConservedExactDiagonalization.jl}
}
```

## Acknowledgments

This package was developed as part of research in quantum many-body physics and exact diagonalization methods.