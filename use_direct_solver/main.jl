# Main program to run AI-REML
# Hongding Gao 2024-03-17
include("xREML.jl")
using .xreml
using Distributions
using Random
using LinearAlgebra
using CSV
using DataFrames
using Statistics
using BenchmarkTools
using SparseArrays



phenodf = CSV.read("pheno.txt", DataFrame; header=false, delim = ' ')
println("Size of phenotype ", size(phenodf))
y = [phenodf[:, 4]; phenodf[:, 5]; phenodf[:, 6]; phenodf[:, 7]; phenodf[:, 8]]
println("Length of y ", length(y))

genodf = CSV.read("SNP.txt", DataFrame; header=false, delim = ' ')
genodf = genodf[:, 2:end]
println("Size of M ", size(genodf))
M0 = Matrix(genodf)
GRMtest = xreml.makeG(M0)
invGRMtest = GRMtest \ I
println("G inverse done")

n = size(phenodf, 1)  
k = size(genodf, 2)  
nt = 5
nind = n

X = sparse([ones(Int8, n) zeros(Int8, n) zeros(Int8, n) zeros(Int8, n) zeros(Int8, n);
            zeros(Int8, n) ones(Int8, n) zeros(Int8, n) zeros(Int8, n) zeros(Int8, n);
            zeros(Int8, n) zeros(Int8, n) ones(Int8, n) zeros(Int8, n) zeros(Int8, n);
            zeros(Int8, n) zeros(Int8, n) zeros(Int8, n) ones(Int8, n) zeros(Int8, n);
            zeros(Int8, n) zeros(Int8, n) zeros(Int8, n) zeros(Int8, n) ones(Int8, n)])

Z = sparse([Matrix{Int8}(1I, n, n) zeros(Int8, n, n) zeros(Int8, n, n) zeros(Int8, n, n) zeros(Int8, n, n);
            zeros(Int8, n, n) Matrix{Int8}(1I, n, n) zeros(Int8, n, n) zeros(Int8, n, n) zeros(Int8, n, n);
            zeros(Int8, n, n) zeros(Int8, n, n) Matrix{Int8}(1I, n, n) zeros(Int8, n, n) zeros(Int8, n, n);
            zeros(Int8, n, n) zeros(Int8, n, n) zeros(Int8, n, n) Matrix{Int8}(1I, n, n) zeros(Int8, n, n);
            zeros(Int8, n, n) zeros(Int8, n, n) zeros(Int8, n, n) zeros(Int8, n, n) Matrix{Int8}(1I, n, n)])

Gprior = Matrix{Float64}(1I, nt, nt)
Rprior = Matrix{Float64}(1I, nt, nt)


@btime xreml.mt_aireml_gblup_new(y, nt, nind, invGRMtest, X, Z, Gprior, Rprior)
@btime xreml.mt_aireml_gblup(y, nt, nind, invGRMtest, X, Z, Gprior, Rprior)









