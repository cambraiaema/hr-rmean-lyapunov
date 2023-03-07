using LightGraphs #importante para matrizes e plotes
using OrdinaryDiffEq
#using DifferentialEquations
using DelimitedFiles
#using Plots; pyplot()

# mudanca uo dentro do for tzise=trt, 0.95

function hr_rede_model!(du, u, p, t)
    
    a, b, c, d, r, s, k, I, N, γ, norm = p[1:11]
    A = p[end]

    du[:, 1] = b.*(u[:, 1].^2) .- a .* (u[:, 1].^3) .+ u[:, 2] .- u[:, 3] .+ I .+ (γ/norm).*(A*u[:, 1])
    du[:, 2] = c .- d.*u[:, 1].^2 .- u[:, 2]
    du[:, 3] = r.*(s.*(u[:, 1] .+ k) .- u[:, 3])

    return nothing
end


function computa_fases(X, Δt, thresh, thresh_B)
    t_size::Int64  = length(X)
    θ::Vector{Float64}       = zeros(Float64, t_size)
    θ_B::Vector{Float64}     = zeros(Float64, t_size)
    tₖ::Vector{Int64} = []
    tₖ_B::Vector{Int64} = []
   # Kₛ = Array{Int64}[]
   # Kₛ_B = Array{Int64}[]
    switch  = true
    switch_B = true
    for t=1:(t_size)
        if (X[t])>thresh
            if switch == false
                #if length(tₖ)>0
                 #   ks = t - tₖ[end]
                 #   Kₛ = vcat(Kₛ, ks)
                #end
                tₖ = vcat(tₖ, t)
            end
            switch = true
        else
            switch = false
        end
#### Bursts
        if (X[t])>thresh_B
            if switch_B == false
               # if length(tₖ_B)>0
               #     ks_B = t - tₖ_B[end]
               #     Kₛ_B = vcat(Kₛ_B, ks_B)
               # end
                tₖ_B = vcat(tₖ_B, t)
            end
            switch_B = true
        else
            switch_B = false
        end
    end
#### Fases
    if length(tₖ)>=2
        for t=1:tₖ[1]
            θ[t] = 2π*(1 + (t - tₖ[1])/(tₖ[2]-tₖ[1]))
        end
        for j=1:(length(tₖ)-1)
            for t=(tₖ[j]+1):tₖ[j+1]
                θ[t] = 2π*(t - tₖ[j])/(tₖ[j+1]-tₖ[j])
            end
        end
        for t=(tₖ[end]+1):t_size
            θ[t] = 2π*(t - tₖ[end])/(tₖ[end]-tₖ[end-1])
        end
    end
    if length(tₖ_B)>=2
        for t=1:tₖ_B[1]
            θ_B[t] = 2π*(1 + (t - tₖ_B[1])/(tₖ_B[2]-tₖ_B[1]))
        end
        for j=1:(length(tₖ_B)-1)
            for t=(tₖ_B[j]+1):tₖ_B[j+1]
                θ_B[t] = 2π*(t - tₖ_B[j])/(tₖ_B[j+1]-tₖ_B[j])
            end
        end
        for t=(tₖ_B[end]+1):t_size
            θ_B[t] = 2π*(t - tₖ_B[end])/(tₖ_B[end]-tₖ_B[end-1])
        end

    end
    return θ, θ_B, tₖ[1], tₖ[end]
end

function resolve_problema(p, u0, tf, tempos)
    
    N           = Int(p[9])
    DS          = ODEProblem(hr_rede_model!, u0, (0.0, tf), p)
    traj  = solve(DS, Tsit5(), dt=0.01, reltol=1e-5, abstol=1e-5, maxiters = 1e7, adaptive=true, saveat=tempos, save_everystep=false)

    X1::Matrix{Float64} = zeros(length(tempos), N)
    
    for t=1:length(tempos)
        X1[t, :] = traj.u[t][(1):(N), 1]
    end

    u0 = traj.u[end]

    return X1, u0

end

function calculaR(t_size, J, Theta, N)
    
    R::Vector{Float64}         = zeros(Float64, t_size*J)

    for idx_T=1:t_size*J
        R[idx_T]=abs(sum(complex.(cos.(Theta[idx_T, 1:N]), sin.(Theta[idx_T, 1:N]))))/N
    end

    return R
end


function main()

    N::Int16                    = 1000
    J::Int16                    = 10
    resolucao::Int16            = 100
    Γ::StepRangeLen             = range(0.0, 0.5, length=resolucao)

    
    a::Float64 = 1.0
    c::Float64 = 1.0
    d::Float64 = 5.0
    r::Float64 = 0.01
    s::Float64 = 4.0
    k::Float64 = 1.6

    norm::Int64    = 12
    A::Matrix{Int8}= readdlm("A_watts_strogatz.txt", '\t', Int, '\n')

    transiente::Float64  = 0
    tf::Float64          = 2000
    Δt::Float64          = 0.5
    tempos::StepRangeLen = range(transiente, stop=tf, step=Δt)
    t_size::Int64        = length(tempos)

    II = [4.39, 4.397, 4.395, 4.29, 4.4, 4.4, 4.4, 4.39, 4.395, 4.395]
    BB = [2.63, 2.625, 2.627, 2.668, 2.623, 2.621, 2.619, 2.627, 2.623, 2.621]

    Threads.@threads for i=1:length(II)

        Rmedio::Vector{Float64}     = zeros(Float64, resolucao)
        Rmedio_B::Vector{Float64}   = zeros(Float64,resolucao)


        b::Float64 = BB[i]
        I::Float64 = II[i]
        
        ctex=rand(-1.5:0.05:1.5, N,1)
        ctey=rand(-12:0.05:0, N,1)
        ctez=rand(3:0.05:4.5, N,1)
        u0::Matrix{Float64}=zeros(N,3)
        u0[:,1]=ctex
        u0[:,2]=ctey
        u0[:,3]=ctez

        for idx_γ=1:resolucao
            
            γ = Γ[idx_γ]
            p = [a, b, c, d, r, s, k, I, N, γ, norm, A]

            x1::Matrix{Float64}  = zeros(Float64, 0,N)
            
            for j=1:J
                X1, u0 = resolve_problema(p, u0, tf, tempos)
                x1 =[x1 ; X1]
            end

            Theta::Matrix{Float64}     = zeros(Float64, t_size*J, N)
            Theta_B::Matrix{Float64}   = zeros(Float64, t_size*J, N)
            Tks_firsts::Vector{Int64}  = zeros(Int64, N)
            Tks_end::Vector{Int64}     = zeros(Int64, N)
            R::Vector{Float64}         = zeros(Float64, t_size*J)
            R_B::Vector{Float64}       = zeros(Float64, t_size*J)

            for idx_N=1:N
                Theta[t_size:end, idx_N], Theta_B[t_size:end, idx_N], Tks_firsts[idx_N], Tks_end[idx_N] = computa_fases(x1[t_size:end,idx_N], Δt, 0.0, 0.88*minimum(x1[t_size:end,idx_N]))
            end
            
            R= calculaR(t_size, J, Theta, N)
            R_B =calculaR(t_size, J, Theta_B, N)
            
            Rmedio[idx_γ] =(1/(minimum(Tks_end)-maximum(Tks_firsts)+1))*sum(R[maximum(Tks_firsts):minimum(Tks_end)])
            Rmedio_B[idx_γ] =(1/(minimum(Tks_end)-maximum(Tks_firsts)+1))*sum(R_B[maximum(Tks_firsts):minimum(Tks_end)])

            open("HR_parord_$(I)_$(b)_smw.txt", "w") do io
                writedlm(io, [Γ, Rmedio_B, Rmedio])
            end

        end
        
    end
end
main()