---
stream: Hamiltonian Scattering Lists
agent: hamiltonian-specialist
started: 2025-09-05T21:00:00Z
status: in_progress
---

## Completed
- Created progress tracking file
- Analyzed existing scattering implementation in main module

## Working On
- Defining dedicated Scattering1Body and Scattering2Body structs
- Creating src/Hamiltonian.jl with scattering framework

## Blocked
- None

## Analysis
The existing implementation in MomentumConservedExactDiagonalization.jl already has:
- Generic Scattering{N} struct (lines 142-156)
- Normal ordering functions for 1-body and 2-body terms (lines 165-217)
- One-body scattering list generation (lines 263-283)
- Two-body scattering list generation with momentum conservation (lines 299-395)
- Integration with Hamiltonian matrix construction (lines 408-460)

However, the requirements call for dedicated Scattering1Body and Scattering2Body structs instead of the generic Scattering{N} approach. The existing implementation is quite comprehensive and already handles:
- Momentum conservation through momentum group pairing
- Upper triangle optimization through normal ordering
- Component conservation
- Hermitian structure maintenance

The main task is to restructure this with the specific struct definitions as requested.