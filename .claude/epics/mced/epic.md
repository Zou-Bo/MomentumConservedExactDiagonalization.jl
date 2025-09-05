---
name: momentum-conserved-ed
status: backlog
created: 2025-09-05T02:41:10Z
updated: 2025-09-05T22:47:08Z
progress: 0%
prd: .claude/prds/mced.md
github: https://github.com/Zou-Bo/MomentumConservedExactDiagonalization.jl/issues/12
---

# Epic: MomentumConservedExactDiagonalization Julia Package

## Overview
Transform existing research code into a production-ready Julia package for momentum-conserved exact diagonalization of quantum many-body systems. The package will provide both core ED capabilities and advanced analysis tools while maintaining momentum conservation symmetries.

### Central Focus
- **Momentum Basis**: Specialized for momentum-space representation of Hilbert space, applicable beyond tight-binding models including Landau levels and moiré systems
- **Flexible Components**: Supports additional components beyond momentum, including both number-conserved and non-conserved components (with hopping terms)
- **Hilbert Space Basis**: Many-Body states are represented with an integer, each bit representing a single-particle basis state
- **Existing Functionality**: All core functions already implemented - focus is on migration, structural organization, and documentation without changing core functionality

## Architecture Decisions

### Core Architecture
- **ED Method Focus**: No real space representations, focus on momentum basis
- **Type System**: Simple MBS <: Integer family for many-body states, no complex structs to save memory usage
- **No Iteration**: ED method has no iteration - direct diagonalization approach
- **Basis Representation**: Hilbert space basis labeled by momentum index and conserved components, mapped to bit positions
- **Hamiltonian Structure**: Scattering lists instead of traditional matrix assembly
- **EDPara-Centric Design**: All parameters, mappings, and scattering amplitudes stored in EDPara container
- **Twisted Boundary Conditions**: Rigid momentum shift support with float (k1, k2) parameters for magnetic flux threading

### Technology Stack
- **Core**: Julia 1.6+ with LinearAlgebra, SparseArrays, KrylovKit
- **Testing**: Test.jl with extensive benchmarks
- **Documentation**: Documenter.jl with physics-focused tutorials
- **CI/CD**: GitHub Actions for cross-platform testing
- **Optional**: Plots.jl for visualization, BenchmarkTools.jl for performance

### Design Patterns
- **EDPara Container**: Single source of truth for all parameters and mappings
- **Scattering Representation**: Hamiltonian terms as lists of scattering processes
- **Integer State Representation**: MBS <: Integer family eliminates need for BasisState/MomentumBasis structs
- **Direct Diagonalization**: No iteration - ED method solves eigenvalue problems directly

## Technical Approach

### Core Components
- **EDPara**: Central parameter container storing momentum indices, conserved components, bit mappings, one-body arrays, scattering amplitude functions, and rigid momentum shift (k1, k2) as Float64 for twisted boundary conditions
- **MBS64 Type Family**: Existing MBS64 <: Integer implementation representing many-body states (up to 64 orbitals) with bit-based occupation
- **Scattering Framework**: Existing scattering-based Hamiltonian construction with 1-body and 2-body terms converted to Scattering structs, with momentum inputs as floats (k_list + momentum_shift)
- **Basis Organization**: Momentum-space Hilbert space construction applicable to Landau levels, moiré systems, and beyond tight-binding models, with support for momentum mesh shifts via rigid translation

### Analysis Modules
- **EntanglementEntropy**: Efficient bipartite entanglement calculations
- **TwistedBoundary**: Magnetic flux threading with spectrum flow tracking
- **TopologicalInvariants**: Many-body Chern numbers via Berry curvature
- **CorrelationAnalysis**: Real-space and momentum-space correlators
- **Visualization**: Basic plotting utilities for spectra and wavefunctions

### Example Systems
- **LandauLevel**: Fractional quantum Hall effect calculations
- **ChernInsulator**: Haldane model implementation
- **MoireSystems**: Simplified twisted bilayer graphene
- **SpinChains**: Heisenberg models with Dzyaloshinskii-Moriya interaction

### Infrastructure
- **Package Structure**: Proper Julia package layout with Project.toml
- **Documentation**: Documenter.jl with physics tutorials
- **Testing**: Comprehensive test suite with known benchmarks
- **CI/CD**: GitHub Actions for Julia 1.6+ compatibility testing

## Implementation Strategy

### Phase 1: Foundation (Weeks 1-3)
- Initialize Julia package structure for existing code migration
- Set up GitHub repository with CI/CD
- Create project skeleton and documentation framework
- Prepare infrastructure for organizing existing EDPara functionality

### Phase 2: Core Engine Organization (Weeks 4-6)
- Migrate existing EDPara parameter container with validation (Task 2)
- Improve EDPara with internal constructor, remove error field, add V_int validation (Task 2a)
- Organize MBS64 type family and bit-based state representation
- Structure existing scattering-based Hamiltonian construction
- Integrate KrylovKit diagonalization routines
- Create comprehensive tests preserving existing functionality

### Phase 3: Analysis Features (Weeks 7-9)
- Implement entanglement entropy calculations
- Add twisted boundary conditions with spectrum flow using rigid momentum shift (k1, k2) in EDPara
- Update interaction functions to accept float momenta (k_list + momentum_shift) instead of integer indices
- Develop many-body Chern number calculations
- Create correlation function utilities

### Phase 4: Polish & Release (Weeks 10-12)
- Complete documentation with tutorials
- Add example systems and physics cases
- Final testing and performance optimization
- Package registration preparation

## Task Breakdown Preview

- [ ] **001.md - Package Foundation**: Initialize Julia package structure and CI/CD (2h, parallel: false)
- [ ] **002.md - Core EDPara Struct**: Migrate existing parameter container with all 9 fields intact (4h, parallel: true)
- [ ] **002a.md - Improved EDPara Struct**: Remove error field, implement internal constructor, add V_int validation (2h, parallel: true)
- [ ] **003.md - Hilbert Space Construction**: Organize existing momentum-resolved basis utilities (6h, parallel: true)
- [ ] **004.md - Hamiltonian Builder**: Structure existing scattering-based Hamiltonian construction (8h, parallel: true)
- [ ] **005.md - Diagonalization Engine**: Integrate existing KrylovKit eigensolvers with momentum conservation (8h, parallel: true)
- [ ] **006.md - Entanglement Entropy**: Implement bipartite entanglement calculations (6h, parallel: true)
- [ ] **007.md - Twisted Boundary Conditions**: Implement magnetic flux threading and spectrum flow (8h, parallel: true)
- [ ] **008.md - Topological Invariants**: Add many-body Chern numbers via Berry curvature (8h, parallel: true)
- [ ] **009.md - Landau Level Example**: Create complete fractional quantum Hall effect example (10h, parallel: true)
- [ ] **010.md - Documentation & Release**: Complete documentation, tutorials, and package registration (12h, parallel: true)

## Tasks Created
- [ ] 001.md - Package Foundation (parallel: false)
- [ ] 002.md - Core EDPara Struct (parallel: true)
- [ ] 002a.md - Improved EDPara Struct (parallel: true)
- [ ] 003.md - Hilbert Space Construction (parallel: true)
- [ ] 004.md - Hamiltonian Builder (parallel: true)
- [ ] 005.md - Diagonalization Engine (parallel: true)
- [ ] 006.md - Entanglement Entropy (parallel: true)
- [ ] 007.md - Twisted Boundary Conditions (parallel: true)
- [ ] 008.md - Topological Invariants (parallel: true)
- [ ] 009.md - Landau Level Example (parallel: true)
- [ ] 010.md - Documentation & Release (parallel: true)

Total tasks: 11
Parallel tasks: 10
Sequential tasks: 1
Estimated total effort: 74 hours (9.25 days)

## Dependencies

### External Services
- **GitHub**: Repository hosting and CI/CD
- **Julia General Registry**: Package distribution
- **Codecov**: Test coverage reporting
- **GitHub Pages**: Documentation hosting

### Internal Dependencies
- **GitHub Repository**: Zou-Bo/MomentumConservedExactDiagonalization.jl
- **Benchmark Data**: Reference results from published papers
- **Testing Infrastructure**: Julia testing ecosystem

### Development Tools
- **Julia**: Core language and package manager
- **Documenter.jl**: Documentation generation
- **BenchmarkTools.jl**: Performance testing
- **JuliaFormatter.jl**: Code formatting

## Success Criteria (Technical)

### Performance Benchmarks
- **Memory**: Handle 16-site systems within 32GB RAM
- **Speed**: 10x faster than naive Python implementations
- **Accuracy**: Energy eigenvalues accurate to 1e-12 relative error
- **Scalability**: Efficient sparse matrix operations for <20 sites

### Quality Gates
- **Test Coverage**: >90% code coverage with meaningful tests
- **Documentation**: Complete API reference with physics background
- **Examples**: 5+ complete physics examples with tutorials
- **Performance**: Regression testing against reference implementations

### Acceptance Criteria
- **Package Registration**: Successfully registered in Julia General
- **User Adoption**: 100+ downloads within 3 months
- **Community**: 3+ external contributors within 6 months
- **Validation**: Results validated against known analytical solutions

## Estimated Effort

### Overall Timeline
- **Total Duration**: 12 weeks (3 months)
- **Development**: 8 weeks core implementation
- **Testing**: 2 weeks comprehensive testing
- **Documentation**: 2 weeks tutorials and examples

### Resource Requirements
- **Primary Developer**: 1 full-time equivalent
- **Testing**: 20% additional time for comprehensive testing
- **Documentation**: 15% additional time for physics tutorials
- **Review**: 5% additional time for code review and optimization

### Critical Path
1. Package foundation and CI/CD setup
2. Core ED engine implementation
3. Advanced analysis features
4. Documentation and examples
5. Package registration and release

## Risk Mitigation

### Technical Risks
- **Memory limitations**: Implement memory-efficient algorithms and clear system size warnings
- **Numerical precision**: Use high-precision arithmetic for degenerate cases
- **Performance bottlenecks**: Profile-driven optimization with BenchmarkTools.jl
- **Julia compatibility**: Test against Julia 1.6 LTS and latest stable

### Development Risks
- **Scope creep**: Strict adherence to PRD scope with feature freeze dates
- **Maintenance burden**: Community transition plan with contributor guidelines
- **User adoption**: Early alpha release for community feedback
- **Documentation quality**: Physics graduate student review and feedback