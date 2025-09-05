---
name: momentum-conserved-ed
status: backlog
created: 2025-09-05T02:41:10Z
progress: 0%
prd: .claude/prds/momentum-conserved-ed.md
github: [Will be updated when synced to GitHub]
---

# Epic: MomentumConservedExactDiagonalization Julia Package

## Overview
Transform existing research code into a production-ready Julia package for momentum-conserved exact diagonalization of quantum many-body systems. The package will provide both core ED capabilities and advanced analysis tools while maintaining momentum conservation symmetries.

## Architecture Decisions

### Core Architecture
- **Modular Design**: Split into core ED engine, analysis modules, and example systems
- **Type System**: Leverage Julia's multiple dispatch for extensible Hamiltonian types
- **Memory Efficiency**: Use sparse matrices and block diagonalization by momentum sectors
- **Performance**: Thread-based parallelism for independent momentum sectors
- **API Design**: Follow Julia conventions with intuitive struct-based interfaces

### Technology Stack
- **Core**: Julia 1.6+ with LinearAlgebra, SparseArrays, Arpack
- **Testing**: Test.jl with extensive benchmarks
- **Documentation**: Documenter.jl with physics-focused tutorials
- **CI/CD**: GitHub Actions for cross-platform testing
- **Optional**: Plots.jl for visualization, BenchmarkTools.jl for performance

### Design Patterns
- **Builder Pattern**: EDPara struct with fluent interface for parameter construction
- **Strategy Pattern**: Pluggable algorithms for different diagonalization methods
- **Template Method**: Common interface for different physical systems
- **Observer Pattern**: Progress tracking for long calculations

## Technical Approach

### Core Components
- **EDPara**: Parameter container with validation and serialization
- **HilbertSpace**: Momentum-resolved basis construction with symmetry exploitation
- **HamiltonianBuilder**: Sparse matrix assembly with block structure
- **Diagonalizer**: Unified interface for different eigensolvers (dense, sparse, iterative)
- **SymmetryAnalyzer**: Automatic symmetry detection and exploitation

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
- Initialize Julia package structure
- Set up GitHub repository with CI/CD
- Create basic project skeleton and documentation framework
- Implement core EDPara struct with validation

### Phase 2: Core Engine (Weeks 4-6)
- Build Hilbert space construction utilities
- Develop Hamiltonian matrix assembly
- Add diagonalization routines with multiple solver options
- Create comprehensive tests and benchmarks

### Phase 3: Analysis Features (Weeks 7-9)
- Implement entanglement entropy calculations
- Add twisted boundary conditions with spectrum flow
- Develop many-body Chern number calculations
- Create correlation function utilities

### Phase 4: Polish & Release (Weeks 10-12)
- Complete documentation with tutorials
- Add example systems and physics cases
- Final testing and performance optimization
- Package registration preparation

## Task Breakdown Preview

- [ ] **001.md - Package Foundation**: Initialize Julia package structure and CI/CD (2h, parallel: false)
- [ ] **002.md - Core EDPara Struct**: Design parameter container with validation (4h, parallel: true)
- [ ] **003.md - Hilbert Space Construction**: Build momentum-resolved basis utilities (6h, parallel: true)
- [ ] **004.md - Hamiltonian Builder**: Implement sparse matrix assembly with block structure (8h, parallel: true)
- [ ] **005.md - Diagonalization Engine**: Add efficient eigensolvers with momentum conservation (8h, parallel: true)
- [ ] **006.md - Entanglement Entropy**: Implement bipartite entanglement calculations (6h, parallel: true)
- [ ] **007.md - Twisted Boundary Conditions**: Implement magnetic flux threading and spectrum flow (8h, parallel: true)
- [ ] **008.md - Topological Invariants**: Add many-body Chern numbers via Berry curvature (8h, parallel: true)
- [ ] **009.md - Landau Level Example**: Create complete fractional quantum Hall effect example (10h, parallel: true)
- [ ] **010.md - Documentation & Release**: Complete documentation, tutorials, and package registration (12h, parallel: true)

## Tasks Created
- [ ] 001.md - Package Foundation (parallel: false)
- [ ] 002.md - Core EDPara Struct (parallel: true)
- [ ] 003.md - Hilbert Space Construction (parallel: true)
- [ ] 004.md - Hamiltonian Builder (parallel: true)
- [ ] 005.md - Diagonalization Engine (parallel: true)
- [ ] 006.md - Entanglement Entropy (parallel: true)
- [ ] 007.md - Twisted Boundary Conditions (parallel: true)
- [ ] 008.md - Topological Invariants (parallel: true)
- [ ] 009.md - Landau Level Example (parallel: true)
- [ ] 010.md - Documentation & Release (parallel: true)

Total tasks: 10
Parallel tasks: 9
Sequential tasks: 1
Estimated total effort: 72 hours (9 days)

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