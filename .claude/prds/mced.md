---
name: mced
description: Julia package for momentum-conserved exact diagonalization of quantum many-body systems with advanced analysis capabilities
status: backlog
created: 2025-09-05T02:41:10Z
---

# PRD: MomentumConservedExactDiagonalization Julia Package

## Executive Summary

This Julia package provides a comprehensive framework for exact diagonalization of quantum many-body systems with momentum conservation, targeting researchers in condensed matter physics. The package transforms existing research code into a production-ready Julia package with advanced features including entanglement entropy calculations, twisted boundary conditions, many-body Chern numbers, and spectrum flow analysis.

### Central Focus
- **Momentum Basis**: Specialized for momentum-space representation of Hilbert space, applicable beyond tight-binding models including Landau levels and moiré systems
- **Flexible Components**: Supports additional components beyond momentum, including both number-conserved and non-conserved components (with hopping terms)
- **Hilbert Space Basis**: Many-Body states are represented with an integer, each bit representing an single-particle basis state.
- **Existing Functionality**: All core functions already implemented - focus is on migration, structural organization, and documentation without changing core functionality

## Problem Statement

Current exact diagonalization implementations for momentum-conserved quantum systems are fragmented, difficult to reproduce, and lack standardized interfaces. Researchers typically maintain custom scripts that are hard to share, validate, or extend. There's a critical need for a unified, well-tested Julia package that provides both core ED capabilities and advanced analysis tools while maintaining momentum conservation symmetries.

## User Stories

### Primary Users
- **Condensed Matter Researchers**: Need to diagonalize Hamiltonians with momentum conservation for studying strongly correlated systems
- **Graduate Students**: Require accessible tools for learning and reproducing ED calculations
- **Computational Physicists**: Need extensible framework for developing new algorithms

### Detailed User Journeys

#### Story 1: Basic ED Calculation
**As a** researcher studying fractional quantum Hall states  
**I want** to diagonalize a Hamiltonian in the lowest Landau level  
**So that** I can obtain energy spectra and eigenstates  

**Acceptance Criteria:**
- Input: Hilbert space definition with 2D momentum labels
- Input: One-body and two-body interaction terms
- Output: Energy eigenvalues and eigenvectors
- Output: Momentum-resolved spectra

#### Story 2: Advanced Analysis
**As a** researcher studying topological properties  
**I want** to calculate entanglement entropy and many-body Chern numbers  
**So that** I can characterize topological phases  

**Acceptance Criteria:**
- Calculate bipartite entanglement entropy for any subsystem
- Compute many-body Chern numbers from twisted boundary conditions
- Generate spectrum flow under boundary condition changes
- Provide visualization tools for results

#### Story 3: Package Usage
**As a** graduate student  
**I want** to install and use the package via Julia's package manager  
**So that** I can reproduce published results and conduct my own research  

**Acceptance Criteria:**
- Simple `] add MomentumConservedExactDiagonalization` installation
- Comprehensive documentation with examples
- Jupyter notebook tutorials
- Benchmarking against known results

## Requirements

### Functional Requirements

#### Core ED Engine (Migration & Structural Organization)
- **Momentum Basis**: Organize existing momentum-space Hilbert space construction, applicable to Landau levels, moiré systems, and beyond tight-binding models
- **EDPara struct**: Centralize existing parameter container with momentum indices, conserved components, bit mappings, one-body arrays, and scattering amplitude functions
- **Component Flexibility**: Structure support for additional conserved quantities including number-conserved and non-conserved components (with hopping terms)
- **MBS64 Type Family**: Organize existing MBS64 <: Integer implementation for many-body states without changing core functionality
- **Scattering Framework**: Structure existing scattering-based Hamiltonian construction with 1-body and 2-body terms
- **Basis Organization**: Organize existing bit-position mapping for momentum and conserved quantum numbers
- **Diagonalization**: Structure existing KrylovKit-based eigensolver without functional changes
- **Threading**: Organize existing threaded Hamiltonian construction while preserving performance

#### Advanced Analysis Features
- **Entanglement entropy**: Calculation for arbitrary bipartitions
- **Twisted boundary conditions**: Implementation of magnetic flux threading
- **Spectrum flow**: Automated tracking of energy levels under boundary changes
- **Many-body Chern numbers**: Calculation via Berry curvature integration
- **Correlation functions**: Real-space and momentum-space correlators
- **Structure factors**: Static and dynamic structure factors

#### Package Infrastructure
- **Julia package structure**: Proper Project.toml, src/, test/, docs/ layout
- **Documentation**: Documenter.jl-based documentation with examples
- **Testing**: Comprehensive test suite with known benchmarks
- **CI/CD**: GitHub Actions for testing across Julia versions
- **Registry submission**: Prepare for General registry inclusion

#### Example Systems
- **Landau level systems**: Fractional quantum Hall effect examples
- **Chern insulators**: Haldane model and variants
- **Moire systems**: Twisted bilayer graphene simplified models
- **Spin systems**: Heisenberg models with Dzyaloshinskii-Moriya interaction

### Non-Functional Requirements

#### Performance
- **Memory efficiency**: Support systems up to ~20 sites on typical hardware
- **Scalability**: Efficient sparse matrix operations using Julia's SuiteSparse
- **Parallelization**: Thread-based parallelization for independent momentum sectors
- **Benchmarking**: Performance regression testing against reference implementations

#### Reliability
- **Numerical stability**: Careful handling of degeneracies and near-degeneracies
- **Validation**: Cross-check results against exact solutions where available
- **Error handling**: Comprehensive error messages for invalid inputs
- **Reproducibility**: Deterministic random number generation for Monte Carlo sampling

#### Usability
- **API design**: Intuitive Julia interface following language conventions
- **Documentation**: Complete docstrings with mathematical background
- **Examples**: Extensive example gallery with physics explanations
- **Tutorials**: Step-by-step guides for common use cases

## Success Criteria

### Technical Metrics
- **Performance**: 10x speedup over naive Python implementations
- **Memory**: Handle 16-site systems within 32GB RAM
- **Accuracy**: Energy eigenvalues accurate to 1e-12 relative error
- **Coverage**: >90% code coverage with meaningful tests

### User Adoption Metrics
- **Downloads**: 100+ downloads within first 3 months
- **Issues**: <5% bug reports vs feature requests
- **Contributions**: 3+ external contributors within 6 months
- **Citations**: 2+ papers citing the package within first year

### Feature Completeness
- **Core features**: All basic ED functionality implemented and tested
- **Advanced features**: 80% of planned analysis tools available
- **Examples**: 5+ complete example systems with tutorials
- **Documentation**: Complete API reference and user guide

## Constraints & Assumptions

### Technical Constraints
- **Julia version**: Support Julia 1.6+ (LTS) and latest stable
- **Dependencies**: Minimal external dependencies beyond LinearAlgebra and SparseArrays
- **Platform**: Cross-platform support (Linux, macOS, Windows)
- **Precision**: Double precision floating point arithmetic
- **Code Structure**: Must follow existing codebase patterns - MBS <: Integer type family, EDPara-centric design, scattering-based Hamiltonian

### Resource Constraints
- **Development time**: 3-month initial development cycle
- **Testing hardware**: Standard academic computing resources
- **Documentation**: Written by physics graduate student target audience
- **Maintenance**: Single maintainer initially, plan for community transition

### Assumptions
- **User expertise**: Graduate-level quantum mechanics knowledge
- **System sizes**: Focus on systems with <20 sites for exact diagonalization
- **Physics domains**: Primarily condensed matter and quantum information
- **Validation**: Users will validate against known analytical results

## Out of Scope

### Explicitly Excluded Features
- **Quantum Monte Carlo**: Stochastic methods beyond exact diagonalization
- **Tensor networks**: MPS/PEPS algorithms for larger systems
- **Real-time dynamics**: Time evolution beyond ground state properties
- **Finite temperature**: Thermodynamic calculations
- **GPU acceleration**: Initial focus on CPU-only implementation
- **Interactive visualization**: Basic plotting only, no GUI
- **Many-body localization**: Specialized algorithms beyond scope
- **Open quantum systems**: Lindblad/master equation approaches

### Future Considerations
- **DMRG integration**: Potential interface with ITensors.jl
- **Machine learning**: Neural network quantum states
- **Experimental comparison**: Direct comparison with experimental data
- **Extended symmetries**: Particle number, spin, and other conservation laws

## Dependencies

### External Dependencies
- **Julia packages**: LinearAlgebra, SparseArrays, KrylovKit
- **Optional packages**: Plots.jl for visualization, Documenter.jl for docs
- **Testing**: Test.jl, Aqua.jl for quality assurance
- **CI/CD**: GitHub Actions, Codecov for coverage reporting

### Internal Dependencies
- **GitHub repository**: Zou-Bo/MomentumConservedExactDiagonalization.jl
- **Documentation hosting**: GitHub Pages for documentation site
- **Registry**: Julia General registry for package distribution
- **Benchmark data**: Reference results from published papers

### Development Dependencies
- **Development tools**: Revise.jl for rapid development
- **Benchmarking**: BenchmarkTools.jl for performance testing
- **Quality assurance**: JET.jl for static analysis, JuliaFormatter.jl
- **Documentation**: Documenter.jl, LiveServer.jl for local docs

## Implementation Timeline

### Phase 1: Package Structure (Week 1-2)
- Initialize Julia package structure
- Set up GitHub repository with CI/CD
- Create basic project skeleton
- Write initial documentation framework

### Phase 2: Core ED Engine Organization (Week 3-6)
- Organize existing EDPara struct with momentum indices and bit mappings
- Structure existing MBS64 <: Integer type family for state representation  
- Organize existing scattering-based Hamiltonian construction
- Structure existing KrylovKit diagonalization routines
- Create comprehensive tests preserving existing functionality

### Phase 3: Advanced Features (Week 7-10)
- Implement entanglement entropy calculations
- Add twisted boundary conditions
- Develop spectrum flow analysis
- Create many-body Chern number calculations
- Write comprehensive example systems

### Phase 4: Documentation & Release (Week 11-12)
- Complete documentation with tutorials
- Final testing and benchmarking
- Package registration preparation
- Community announcement and outreach

## Risk Assessment

### Technical Risks
- **Memory limitations**: May need to limit system sizes
- **Numerical precision**: Degenerate eigenspaces could cause issues
- **Performance bottlenecks**: Sparse matrix operations might be slow
- **Julia compatibility**: Changes in Julia 1.x might break code

### Mitigation Strategies
- **Incremental development**: Test with small systems first
- **Extensive testing**: Benchmark against known results
- **Community feedback**: Early alpha release for testing
- **Documentation**: Clear limitations and expected performance

### Success Risks
- **User adoption**: Package might be too specialized
- **Maintenance burden**: Single maintainer model unsustainable
- **Competition**: Existing packages might already solve the problem
- **Scope creep**: Adding too many features too quickly