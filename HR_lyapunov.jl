using LaTeXStrings
using LightGraphs
using DynamicalSystems
using DelimitedFiles
#
@inline @inbounds function hr_model(u, p, t)
    
    a, b, c, d, r, s, k, I = p
    z1, z2, z3 = u

    dz1 = b*(z1^2)-a*(z1^3)+z2-z3+I
    dz2 = c-d*z1^2-z2
    dz3 = r*(s*(z1+k)-z3)

    return SVector{3}(dz1, dz2, dz3)
end

function main()
    resolucao = 1000

    Δt          = 0.01
    u0          = [0.1, 0.2, 0.3]
    Tt          = 5000

    Ii          = range(4.25, 4.45, length=resolucao)
    B           = range(2.6, 2.7, length=resolucao)
    Exp         = zeros(resolucao, resolucao)
#
    Threads.@threads for idx_I=1:resolucao
        for idx_b=1:resolucao
            a = 1.0
            b = B[idx_b]
            c = 1.0
            d = 5.0
            r = 0.01
            s = 4.0
            k = 1.6
            I = Ii[idx_I]

            p                   = [a; b; c; d; r; s; k; I]
            prob                = ContinuousDynamicalSystem(hr_model, u0, p)
            Exp[idx_I, idx_b]   = lyapunovspectrum(prob, Tt, Δt = 0.01; Ttr=25000.0)[1]
        end

        open("HR_lya5kTtr25k.txt", "w") do io
            writedlm(io, [B Ii Exp])
        end
        # println(idx_I/resolucao*100)
    end
    nothing
end

main()