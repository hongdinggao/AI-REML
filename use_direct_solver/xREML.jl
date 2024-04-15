# functions for the augmented AI-REML (mt_aireml_gblup_new) and standard AI-REML (mt_aireml_gblup)
# based on multi-trait GBLUP model using direct solving method
# Hongding Gao 2024-03-17
module xreml
  using Distributions
  using Random
  using LinearAlgebra
  using CSV
  using DataFrames
  using Statistics
  using SparseArrays


function mul_offdiag2(A::Matrix)
    n = size(A, 1)
    for i in 1:n
        for j in (i+1):n  
            A[i, j] *= 2
            A[j, i] *= 2  
        end
    end
    return A
end


function mat2vec(M::AbstractMatrix{T}) where T
    n = size(M, 1)
    v = Vector{T}(undef, n*(n+1)÷2)  
    k = 1
    for i in 1:n
        for j in i:n
            v[k] = M[i, j]
            k += 1
        end
    end
    return v
end


function vec2mat(v::Vector{T}) where T
    n = length(v)
    m = Int(sqrt(2n + 0.25) - 0.5)
    M = Matrix{T}(undef, m, m)
    k = 1
    for i in 1:m
        for j in i:m
            M[i, j] = v[k]
            M[j, i] = v[k]  
            k += 1
        end
    end
    return M
end
  
# AI REML for multi-trait GBLUP model
# GRM is the genomic relationship matrix while G0 is the genetic VCV matrix
function mt_aireml_gblup(y, ntrait, nind, invGRM, X, Z, G0, R0, tol=1.0e-12)
    n = div(length(y), ntrait)             
    rankX = ntrait                        # assume only general mean as the fixed effects
    rankA = ntrait * nind
    W = [X Z]
    nmk = nind       
    nfix = 1         
    u = zeros(nmk, ntrait)           
    e = zeros(n*ntrait)
    sol = zeros(rankX+rankA)         
    neq = length(sol)                
    z = rankX + 1     
    oldG = copy(G0)
    oldR = copy(R0)
    G = Array{Float64}(undef, ntrait, ntrait)
    R = Array{Float64}(undef, ntrait, ntrait)
    gG = Array{Float64}(undef, ntrait, ntrait)
    rG = Array{Float64}(undef, ntrait, ntrait)
    diff = 1
    i = 0
    while diff > tol
        i += 1
        iG0 = inv(oldG)
        iR0 = inv(oldR)
        Ginv = kron(iG0, invGRM)
        Rinv = sparse(kron(iR0, Matrix(1.0I, n, n)))
        tX_Rinv = X'Rinv
        tZ_Rinv = Z'Rinv
        tW_Rinv = [tX_Rinv;
                   tZ_Rinv]
        RHS = tW_Rinv * y
        LHS = Matrix([tX_Rinv * X tX_Rinv * Z;
                      tZ_Rinv * X tZ_Rinv * Z .+ Ginv])
        invLHS = LHS \ I
        sol .= invLHS * RHS
        e .= y .- W * sol
        u .= reshape(sol[z:end], nmk, ntrait)
        f2G = u'*invGRM*u
        e2 = reshape(e, n, ntrait)
        f2R = e2'e2
        f1G = zeros(ntrait, ntrait)
        c22 = invLHS[z:end, z:end]
        T = kron(Matrix(1.0I, ntrait, ntrait), invGRM) * c22
        for j ∈ 1:nmk
            idx = j .+ [0:(ntrait-1);] .* nmk 
            f1G .= f1G .+ T[idx, idx]
        end
        f1R = zeros(ntrait, ntrait)
        PEV = W * invLHS * W'
        for k ∈ 1:n
            idx = k .+ [0:(ntrait-1);] .* n
            f1R .= f1R .+ PEV[idx, idx]
        end

        bigU = u * iG0
        Rinv_e = Rinv * e
        npar = ntrait * (ntrait+1)
        F = zeros(n*ntrait, npar)
        T = zeros(neq, npar)
        ij = 1
        smallnpar = npar ÷ 2
        for j ∈ 1:ntrait
            jjo = [1:n;] .+ (j-1)*n
            jja = [1:nmk;] .+ (j-1)*nmk
            for i ∈ j:ntrait
                iio = [1:n;] .+ (i-1)*n
                iia = [1:nmk;] .+ (i-1)*nmk
                if (i == j)
                    F[:, ij] .= Z[:, iia] * bigU[:, i]
                    RHS = tW_Rinv * F[:, ij]
                    T[:, ij] .= invLHS * RHS
                    D = zeros(Int8, n*ntrait, n*ntrait)
                    D[iio, iio] = Matrix{Int8}(1I, n, n)
                    F[:, ij+smallnpar] .= D * Rinv_e
                    RHS = tW_Rinv * F[:, ij+smallnpar]
                    T[:, ij+smallnpar] .= invLHS * RHS
                else
                    F[:, ij] .= Z[:, iia] * bigU[:, j] .+ Z[:, jja] * bigU[:, i]
                    RHS = tW_Rinv * F[:, ij]
                    T[:, ij] .= invLHS * RHS
                    D = zeros(Int8, n*ntrait, n*ntrait)
                    D[iio, jjo] = Matrix{Int8}(1I, n, n)
                    D[jjo, iio] = Matrix{Int8}(1I, n, n)
                    F[:, ij+smallnpar] .= D * Rinv_e
                    RHS = tW_Rinv * F[:, ij+smallnpar]
                    T[:, ij+smallnpar] .= invLHS * RHS
                end
                ij += 1 
            end
        end
        
        Iai = inv(F'*Rinv*F .- T'*tW_Rinv*F)
        gG = iG0 * (nmk .* oldG .- f1G .- f2G) * iG0
        gG = mul_offdiag2(gG)
        rG = iR0 * (n .* oldR .- f1R .- f2R) * iR0
        rG = mul_offdiag2(rG)
        Grad = [mat2vec(gG); mat2vec(rG)]
        thetaAI = [mat2vec(oldG); mat2vec(oldR)] .- Iai*Grad
        theta = thetaAI

        G = vec2mat(theta[begin:(ntrait*(ntrait+1)÷2)])
        R = vec2mat(theta[(ntrait*(ntrait+1)÷2+1):end])

        println(i)
        println(G)
        println(R)
        oldtheta = [mat2vec(oldG); mat2vec(oldR)]
        diff = ([theta .- oldtheta]' * [theta .- oldtheta])/(theta'theta)
        oldG = copy(G)
        oldR = copy(R)
    end
    println(i)

    return (G, R)
end # of mt_aireml_gblup

# New AI REML for multi-trait GBLUP model
# GRM is the genomic relationship matrix while G0 is the genetic VCV matrix
function mt_aireml_gblup_new(y, ntrait, nind, invGRM, X, Z, G0, R0, tol=1.0e-12)
    n = div(length(y), ntrait)             
    rankX = ntrait
    rankA = ntrait * nind
    W = [X Z]
    N_EQ = rankX + rankA
    nmk = nind
    nfix = 1         
    u = zeros(nmk, ntrait)           
    bigU = similar(u)
    e = zeros(n*ntrait)
    sol = zeros(rankX+rankA)         
    z = rankX + 1     
    oldG = copy(G0)
    oldR = copy(R0)
    G = Array{Float64}(undef, ntrait, ntrait)
    R = Array{Float64}(undef, ntrait, ntrait)
    gG = Array{Float64}(undef, ntrait, ntrait)
    rG = Array{Float64}(undef, ntrait, ntrait)
    diff = 1
    i = 0
    while diff > tol 
        i += 1
        iG0 = inv(oldG)
        iR0 = inv(oldR)
        Ginv = kron(iG0, invGRM)
        Rinv = sparse(kron(iR0, Matrix(1.0I, n, n)))
        tX_Rinv = X'Rinv
        tZ_Rinv = Z'Rinv
        tW_Rinv = [tX_Rinv; 
                   tZ_Rinv]
        iRy = Rinv * y
        RHS = tW_Rinv * y
        LHS = Matrix([tX_Rinv * X tX_Rinv * Z;
                      tZ_Rinv * X tZ_Rinv * Z .+ Ginv])
        invLHS = LHS \ I
        sol .= invLHS * RHS
        e .= y .- W * sol
        u .= reshape(sol[z:end], nmk, ntrait)
        f1G = zeros(ntrait, ntrait)
        c22 = invLHS[z:end, z:end]
        T = kron(Matrix(1.0I, ntrait, ntrait), invGRM) * c22 
        for j ∈ 1:nmk
            idx = j .+ [0:(ntrait-1);] .* nmk 
            f1G .= f1G .+ T[idx, idx]
        end
        f1R = zeros(ntrait, ntrait)
        PEV = W * invLHS * W'
        for k ∈ 1:n
            idx = k .+ [0:(ntrait-1);] .* n
            f1R .= f1R .+ PEV[idx, idx]
        end
        
        Rinv_e = Rinv * e
        mul!(bigU, u, iG0)
        npar = ntrait * (ntrait+1)
        F = zeros(n*ntrait, npar)
        ij = 1
        smallnpar = npar ÷ 2
        for j ∈ 1:ntrait
            jjo = [1:n;] .+ (j-1)*n
            jja = [1:nmk;] .+ (j-1)*nmk
            for i ∈ j:ntrait
                iio = [1:n;] .+ (i-1)*n
                iia = [1:nmk;] .+ (i-1)*nmk
                if (i == j)
                    F[iio, ij] .= Z[iio, iia] * bigU[:, i]
                    F[iio, ij+smallnpar] .= Rinv_e[iio]
                else
                    F[jjo, ij] .= Z[jjo, jja] * bigU[:, i]
                    F[iio, ij] .= Z[iio, iia] * bigU[:, j]
                    F[jjo, ij+smallnpar] .= Rinv_e[iio]
                    F[iio, ij+smallnpar] .= Rinv_e[jjo]
                end
                ij += 1 
            end
        end

        gG = iG0 * (nmk .* oldG .- f1G) * iG0
        gG = mul_offdiag2(gG)
        rG = iR0 * (n .* oldR .- f1R) * iR0
        rG = mul_offdiag2(rG)
        Grad = [mat2vec(gG); mat2vec(rG)]
        fMf = zeros(npar, npar)
        fMf .= F' * Rinv * F .- F' * Rinv * PEV * Rinv * F
        fMy = zeros(npar)
        fMy .= F' * iRy .- Grad - F' * Rinv * PEV * iRy 
        sol_new = fMf \ fMy
        
        thetaAI = [mat2vec(oldG); mat2vec(oldR)] .+ sol_new
        theta = thetaAI
        
        G = vec2mat(theta[begin:(ntrait*(ntrait+1)÷2)])
        R = vec2mat(theta[(ntrait*(ntrait+1)÷2+1):end])
        println(i)
        println(G)
        println(R)
        oldtheta = [mat2vec(oldG); mat2vec(oldR)]
        diff = ([theta .- oldtheta]' * [theta .- oldtheta])/(theta'theta)
        oldG = copy(G)
        oldR = copy(R)
    end
    println(i)
    
    return (G, R)
end # of mt_aireml_gblup_new


function makeG(M)
    n = size(M, 1)
    G = Array{Float64}(undef, n, n)
    p = mean(M, dims=1)/2
    q = 1 .- p
    pq = [p; q]
    maf = minimum(pq, dims=1)
    maf = maf .* 2
    M = M .- maf
    sigma2pq = sum(2*p.*q)
    G = M*M'/sigma2pq
    return(G)
end



end # of the module
