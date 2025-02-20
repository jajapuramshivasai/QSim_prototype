using BenchmarkTools
# using .QSim_opt  # Ensure your module is properly imported
using SparseArrays, LinearAlgebra

function create_ghz(n::Int)
    # Initialize n-qubit state |0⟩^⊗n
    s = QSim_opt.sv(n, 0)
    
    # Create superposition
    QSim_opt.h!(s, 1)
    
    # Create entanglement
    for target in 2:n
        QSim_opt.cnot!(s, 1, target)
    end
    
    s
end

function benchmark_ghz(n_max::Int, step::Int=1)
    results = []
    for n in 1:step:n_max
        # Benchmark state preparation
        bench = @benchmarkable create_ghz($n) samples=5 evals=1
        t = run(bench)
        
        # Verify state correctness
        s = create_ghz(n)
        valid = true
        expected = 1/√2
        for i in 0:(1<<n)-1
            amp = s[i+1]
            if (i == 0 || i == (1<<n)-1)
                valid &= abs(amp - expected) < 1e-4
            else
                valid &= abs(amp) < 1e-4
            end
        end
        
        # Memory statistics
        mem_alloc = t.memory
        mem_usage = sizeof(s)
        
        push!(results, (n=n, time=median(t.times)/1e9, 
                       memory=mem_alloc/1e6, valid=valid))
        
        println("n=$n: $(round(results[end].time, digits=3)) s, ",
                "Memory: $(round(mem_usage/1e6, digits=2)) MB, ",
                "Valid: $valid")
    end
    results
end

# Run benchmark up to 25 qubits (adjust based on your system)
benchmark_results = benchmark_ghz(25, 5)

# Optional: Plotting (requires Plots.jl)
# using Plots
# plot([r.n for r in results], [r.time for r in results], 
#      xlabel="Number of Qubits", ylabel="Time (s)", 
#      title="GHZ State Preparation Time", ma
