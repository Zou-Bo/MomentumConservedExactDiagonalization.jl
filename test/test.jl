using Test
using MomentumConservedExactDiagonalization

k_list = [0 1 2 0 1 2 0 1 2;
          0 0 0 1 1 1 2 2 2]
Nk = 9; Gk = (3,3)
m = 3; @assert iszero(mod(Nk, m))
Ne = Nk ÷ m
Gl = sqrt(2π/sqrt(0.75))
D_l = 5.0
W0 = 1.0
function VFF(q1::Float64, q2::Float64)
    ql = sqrt(q1^2 + q2^2 - q1*q2) * Gl
    if ql == 0.0
        return W0 * D_l 
    end
    return W0 / ql * tanh(ql * D_l) * exp(-0.5 * ql^2)
end
function ita(g1::Int64, g2::Int64)
    if iseven(g1) && iseven(g2)
        return 1
    else
        return -1
    end
end
function ql_cross(q1_1, q1_2, q2_1, q2_2)
    return q1_1 * q2_2 - q1_2 * q2_1
end
function V_int(kf1, kf2, ki1, ki2, cf1=1, cf2=1, ci1=1, ci2=1; output=false)::ComplexF64
    ki1 = k_list[1:2, ki1]
    ki2 = k_list[1:2, ki2]
    kf1 = k_list[1:2, kf1]
    kf2 = k_list[1:2, kf2]
    q = rem.(ki1 - kf1, Gk, RoundNearest)
    G_shift1 = (ki1 - kf1 - q) .÷ Gk
    G_shift2 = (kf2 - ki2 - q) .÷ Gk

    ki1 = ki1 ./ Gk
    ki2 = ki2 ./ Gk
    kf1 = kf1 ./ Gk
    kf2 = kf2 ./ Gk

    V_total = 0.0 + 0.0im
    # N shells of reciprocal lattice vectors G
    Nshell = 2
    for g1 in -Nshell:Nshell, g2 in -Nshell:Nshell
        if abs(g1-g2) > Nshell
            continue
        end

        qq1 = q[1] / Gk[1] + g1
        qq2 = q[2] / Gk[2] + g2

        phase_angle = 0.5ql_cross(ki1[1], ki1[2], kf1[1], kf1[2])
        phase_angle += 0.5ql_cross(ki1[1]+kf1[1], ki1[2]+kf1[2], qq1, qq2)
        phase_angle += 0.5ql_cross(ki2[1], ki2[2], kf2[1], kf2[2])
        phase_angle += 0.5ql_cross(ki2[1]+kf2[1], ki2[2]+kf2[2], -qq1, -qq2)

        phase = cispi(2.0phase_angle)
        sign = ita(g1+G_shift1[1], g2+G_shift1[2]) * ita(g1+G_shift2[1], g2+G_shift2[2])

        V_total += sign * phase * VFF(qq1, qq2)
    end

    return V_total
end
para = EDPara(k_list=k_list, Gk=Gk, V_int = V_int);
blocks, block_k1, block_k2, k0number = 
    ED_momentum_block_division(para, ED_mbslist(para, (Ne,)));
length.(blocks)
scat_list1 = ED_sortedScatteringList_onebody(para);
scat_list2 = ED_sortedScatteringList_twobody(para);
@time EDsolve(blocks[1], scat_list1, scat_list2, 10)


@test 1==1