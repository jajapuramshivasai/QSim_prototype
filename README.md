# QSim Prototype

## Overview
QSim Prototype is a high-performance quantum simulator designed for hybrid tasks utilizing parallel computing. The project aims to provide an efficient and scalable framework for simulating quantum circuits, with a focus on performance optimization and parallel computation.

## Features
- **High-Performance Quantum Simulation**: Optimized for speed and efficiency, leveraging parallel computing techniques.
- **Modular Design**: The simulator is structured into modules, each handling specific aspects of quantum computation.
- **Gate Operations**: Supports a variety of quantum gates, including Hadamard, Pauli-X, Pauli-Y, Pauli-Z, and controlled gates like CNOT.
- **State Vector Manipulation**: Efficient state vector operations for simulating quantum states.
- **Benchmarking Tools**: Includes tools for benchmarking the performance of quantum circuits.
- **Custom Circuit Creation**: Allows users to define and simulate custom quantum circuits.

## Directory Structure
  ├── README.md
  
  ├── LICENSE
  
  ├── QSim_opt
  
  ├── benchmark.txt
  
  ├── circuit.jl
  
  ├── code_native
  
  ├── fin.jl
  
  ├── gates.jl
  
  ├── kronnecker_prod.jl
  
  ├── temp2.jl
  
  └── to_do.txt

## Modules
- **QSim_opt**: Core module containing optimized quantum operations and state vector manipulations.
- **circuit.jl**: Implements quantum circuit creation and manipulation.
- **fin.jl**: Contains the main logic for quantum circuit simulation, including state vector initialization and gate application.
- **gates.jl**: Defines various quantum gates and their operations.
- **kronnecker_prod.jl**: Implements Kronecker product operations for tensor networks. (Currently Under Progress)
- **temp2.jl**: Benchmarking tools for testing the performance of quantum circuits.

## Usage
To use the QSim Prototype, you can start by defining a quantum circuit and applying gates to it. Here's an example of creating a GHZ state:

```julia
using .QSim_opt

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

# Create a 5-qubit GHZ state
ghz_state = create_ghz(5)
```

## Benchmarking
The project includes benchmarking tools to measure the performance of quantum circuits. You can run benchmarks to test the time and memory usage of different quantum operations.

```julia
using BenchmarkTools

function benchmark_ghz(n_max::Int, step::Int=1)
    results = []
    for n in 1:step:n_max
        bench = @benchmarkable create_ghz($n) samples=5 evals=1
        t = run(bench)
        push!(results, (n=n, time=median(t.times)/1e9))
        println("n=$n: $(round(results[end].time, digits=3)) s")
    end
    results
end

# Run benchmark up to 10 qubits
benchmark_results = benchmark_ghz(10, 2)
```

## License
This project is licensed under the Apache License 2.0. See the [LICENSE](LICENSE) file for more details.
