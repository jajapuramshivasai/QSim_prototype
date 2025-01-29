#problem: time isnt improving seems like multithreading is not working on mac
#issue:replace "I" with "-" #status:fixed
#issue:replace "⨁" with "⊕" or "⨁" #status:fixed
#issue:change fp64 to fp32 or fp16
#issue:replace c = A*B with mul!(C, A, B) #
#issue:add @inbounds to all for loops
#task: implemnt universal instruction set {U,CU,MCU} and if U coressponds to any known gate then print them as known gate in print_circuit_new
#>Example: U

#implement HQNN and quantum transformer
using BenchmarkTools
using Profile
using LinearAlgebra
using Plots
# using DistributedArrays
using SparseArrays
using Hwloc
# Hwloc.num_physical_cores()
# BLAS.set_num_threads(4)
"""
Hwloc.num_physical_cores() = 8
nthreads() = 1
pennylane limit 26 qubits , Time = 20 seconds
julia naive limit 14 qubits, Time = 38 seconds


naive approach 

qubit_count | complexity(single core)                                          | complexity(parallel 4 cores)+optimized
---------------------------------------------------------------------------------------------------------------------------------------------------
6           | 0.000426 seconds (266 allocations: 804.758 KiB)                  |     
7           | 0.000704 seconds (343 allocations: 3.437 MiB)                    |         
8           | 0.002305 seconds (419 allocations: 15.020 MiB)                   |
9           | 0.013136 seconds (500 allocations: 65.312 MiB, 25.68% gc time)   |
10          | 0.143950 seconds (586 allocations: 282.401 MiB, 71.70% gc time)  |
11          | 0.591260 seconds (677 allocations: 1.186 GiB, 67.93% gc time)    |
12          | 1.744083 seconds (773 allocations: 5.077 GiB, 44.98% gc time)    |
13          | 9.349812 seconds (874 allocations: 21.641 GiB, 13.41% gc time)   |
14          | 40.137828 seconds (980 allocations: 91.892 GiB, 4.96% gc time)


"""

struct QuantumCircuit
    gates::Vector{Vector{String}}
    num_qubits::Int
    function QuantumCircuit(num_qubits::Int)
        circuit = new(Vector{Vector{String}}(), num_qubits)
        add_layer!(circuit)
        return circuit
    end
end

using CSV
using DataFrames

function save_circuit_to_csv(circuit::QuantumCircuit, filename::String)
    data = []
    for (layer_idx, layer) in enumerate(circuit.gates)
        for (qubit_idx, gate) in enumerate(layer)
            push!(data, (layer_idx, qubit_idx, gate))
        end
    end
    df = DataFrame(data, [:Layer, :Qubit, :Gate])
    CSV.write(filename, df)
end
function load_circuit_from_csv(filename::String)::QuantumCircuit
    df = CSV.read(filename, DataFrame)
    num_qubits = maximum(df.Qubit) + 1
    circuit = QuantumCircuit(num_qubits)
    for row in eachrow(df)
        layer = row.Layer
        qubit = row.Qubit
        gate = row.Gate
        set_gate!(circuit, layer, qubit, gate)
    end
    return circuit
end

function add_layer!(circuit::QuantumCircuit)
    push!(circuit.gates, fill("-", circuit.num_qubits)) #fix 1
end

function set_gate!(circuit::QuantumCircuit, layer::Int, qubit::Int, gate::String)
    if layer > length(circuit.gates)
        for _ in 1:(layer - length(circuit.gates))
            add_layer!(circuit)
        end
    end
    circuit.gates[layer][qubit] = gate
end

function H(circuit::QuantumCircuit,qubit::Int)
    layer = length(circuit.gates)
 
    Condition::Bool = (circuit.gates[layer][qubit] == "-") #fix 2
    if Condition
        circuit.gates[layer][qubit] = "H"
    else
        add_layer!(circuit)
        circuit.gates[layer+1][qubit] = "H"
    end

end

function X(circuit::QuantumCircuit,qubit::Int)
    layer = length(circuit.gates)
 
    Condition::Bool = (circuit.gates[layer][qubit] == "-") #fix 2
    if Condition
        circuit.gates[layer][qubit] = "X"
    else
        add_layer!(circuit)
        circuit.gates[layer+1][qubit] = "X"
    end

end

function R(circuit::QuantumCircuit,qubit::Int, theta::Float64, phi::Float64, lambda::Float64)

    R :: Matrix{Complex{Float64}} = [cos(theta/2) -im*exp(im*lambda)*sin(theta/2); im*exp(im*phi)*sin(theta/2) exp(im*(phi+lambda))*cos(theta/2)]
    layer = length(circuit.gates)
 
    Condition::Bool = (circuit.gates[layer][qubit] == "-") #fix 3
    if Condition
        circuit.gates[layer][qubit] = "R{$theta,$phi,$lambda}"
    else
        add_layer!(circuit)
        circuit.gates[layer+1][qubit] = "R($theta,$phi,$lambda)"
    end

end


function CNOT(circuit::QuantumCircuit,control::Int, target::Int)
    layer = length(circuit.gates)
    d = abs(control - target)
    if control < target
        to_apply_nxt_layer = false
        for i in control:target
            if circuit.gates[layer][i] != "-" #fix 4
                to_apply_nxt_layer = true
                break
            end
        end
        if to_apply_nxt_layer
            add_layer!(circuit)
            circuit.gates[layer+1][target] = "⊕$d"
            circuit.gates[layer+1][control] = "●$d"
            for i in control+1:target-1
                circuit.gates[layer+1][i] = "|"
            end
        else
            circuit.gates[layer][target] = "⊕$d"
            circuit.gates[layer][control] = "●$d"
            for i in control+1:target-1
                circuit.gates[layer][i] = "|"
            end
        end

    elseif control > target
        to_apply_nxt_layer = false
        for i in target:control
            if circuit.gates[layer][i] != "-" #fix 5
                to_apply_nxt_layer = true
                break
            end
        end
        if to_apply_nxt_layer
            add_layer!(circuit)
            circuit.gates[layer+1][control] = "●$d"
            circuit.gates[layer+1][target] = "⊕$d"
            for i in target+1:control-1
                circuit.gates[layer+1][i] = "|"
            end
        else
            circuit.gates[layer][control] = "●$d"
            circuit.gates[layer][target] = "⊕$d"
            for i in target+1:control-1
                circuit.gates[layer][i] = "|"
            end
        end

    else
        # println("Control and target qubits are same")
        error("Control and target qubits are same")
    end
end



function identity_(qubits::Int)::Matrix{Complex{Float64}}
    n = 1 << qubits
    return Matrix{Complex{Float64}}(I, n, n)
end

function CX(d::Int)::Matrix{Complex{Float64}} #correst >verified
    P0::Matrix{Complex{Float64}} = [1 0; 0 0]
    P1::Matrix{Complex{Float64}} = [0 0; 0 1]
    X::Matrix{Complex{Float64}} = [0 1; 1 0]
 
    U = kron(P0 , identity_(d) ) + kron( P1 , kron(identity_(d-1) ,X))
    return U
    
end
function XC(d::Int)::Matrix{Complex{Float64}} #correst >verified
    P0::Matrix{Complex{Float64}} = [1 0; 0 0]
    P1::Matrix{Complex{Float64}} = [0 0; 0 1]
    X::Matrix{Complex{Float64}} = [0 1; 1 0]
    U =   kron(identity_(d) ,P0) + kron(X , kron(identity_(d-1) , P1))
    return U
    
end
#new proposed instruction set

function print_circuit(circuit::QuantumCircuit)
    max_qubit_digits = length(string(circuit.num_qubits))
    @inbounds for qubit in 1:circuit.num_qubits
        qubit_str = lpad("Q$qubit:", max_qubit_digits + 2)
        print(qubit_str, " ")
        @inbounds for layer in circuit.gates
            print(layer[qubit][1], "--")
        end
        println()
    end
end

function print_matrix(matrix::Matrix{Complex{Float64}})
    for row in 1:size(matrix, 1)
        for col in 1:size(matrix, 2)
            real_part = real(matrix[row, col])
            imag_part = imag(matrix[row, col])
            if imag_part >= 0
                print("$(real_part)+$(imag_part)i ")
            else
                print("$(real_part)$(imag_part)i ")
            end
        end
        println()
    end
end

function print_state(state_vector::Vector{ComplexF64})
    n = Int(log2(length(state_vector)))
    println("State representation:")
    @inbounds for i in 0:(length(state_vector)-1)
        amplitude = state_vector[i+1]
        if abs(amplitude) > 1e-10  # Only print non-zero amplitudes
            binary = lpad(string(i, base=2), n, '0')
            reversed_binary = reverse(binary)
            println("|$reversed_binary⟩: $amplitude")
        end
    end
end

# function plot_probabilities(state_vector::Vector{ComplexF64})
#     n = Int(log2(length(state_vector)))
#     probabilities = abs2.(state_vector)
#     labels = [reverse(lpad(string(i, base=2), n, '0')) for i in 0:(length(state_vector)-1)]
#     bar(labels, probabilities, xlabel="State", ylabel="Probability", title="State Probabilities")
# end
function plot_probabilities(state_vector::Vector{ComplexF64},name::String="State Probabilities")
    n = Int(log2(length(state_vector)))
    println("State representation:")
    probabilities = []
    labels = []
    @inbounds for i in 0:(length(state_vector)-1)
        amplitude = state_vector[i+1]
        if abs(amplitude) > 1e-10  # Only print non-zero amplitudes
            binary = lpad(string(i, base=2), n, '0')
            reversed_binary = reverse(binary)
            # println("|$reversed_binary⟩: $amplitude")
            # probabilities.push(abs(amplitude)^2)
            push!(probabilities, abs(amplitude)^2)
            push!(labels, "|$reversed_binary⟩")
            # labels.push("|$reversed_binary⟩")
        end
    end
    bar(labels, probabilities, xlabel="State", ylabel="Probability", title=name, legend=false)
end

#make a distributed version of compute_layer

function compute_layer(circuit::QuantumCircuit, layer::Int) #test and verify this function >> Verified
    U::Matrix{Complex{Float64}} = I(1)   #add R gate to this
    qubit = 1
    @inbounds while qubit <= circuit.num_qubits
        gate = circuit.gates[layer][qubit]
        if gate == "-"
            U = kron(U, [1 0; 0 1])
            
        elseif gate == "X"
            U = kron(U, [0 1; 1 0])
        elseif gate == "H"#fix 6
            U = kron(U, [1 1; 1 -1]/sqrt(2))

        elseif startswith(gate, "R(") #Verify this
            params = split(gate[3:end-1], ",")
            theta = parse(Float64, params[1])
            phi = parse(Float64, params[2])
            lambda = parse(Float64, params[3])
            R = [cos(theta/2) -im*exp(im*lambda)*sin(theta/2); im*exp(im*phi)*sin(theta/2) exp(im*(phi+lambda))*cos(theta/2)]
            U = kron(U, R)
        
        elseif gate[1] == '⨁'
            g = parse(Int, gate[4:end]) # change 2 to 4 to make it work
            d = g
            U = kron(U, XC(d))
            qubit += d
        elseif gate[1] == '●'
            g = parse(Int, gate[4:end]) 
            d = g
            U = kron(U, CX(d))
            qubit += d
        else
            error("Unknown gate: $gate")
        end
        qubit += 1
    end
    return U
end



function initialize_statevector(n::Int, m::Int)::Vector{ComplexF64}
    size = 1<<n
    result = zeros(ComplexF64, size)
    result[m+1] = 1.0
    return result
end
function compute_sv(circuit::QuantumCircuit, state::Vector{ComplexF64})
    @inbounds for layer in 1:length(circuit.gates)
        U = compute_layer(circuit, layer)
        state .= U * state
        # mul!(state,U,state)
    end
    
end

function compute_circ(circuit::QuantumCircuit)
    A::Matrix{ComplexF64} = I(1<<circuit.num_qubits)
    @inbounds for layer in 1:length(circuit.gates)
        U = compute_layer(circuit, layer)
        A .= U * A
        # mul!(A,U,A)
    end
    return A
    
end
function get_memory_size(obj)
    # Use @allocated macro to measure memory allocation during evaluation
    mem = @allocated begin
        obj
    end
    return mem
end

# Example usage
# num_qubits =2
# state = initialize_statevector(num_qubits, 0)
# qc = QuantumCircuit(num_qubits)
# set_gate!(qc, 1, 1, "H")
# set_gate!(qc, 1, 2, "-")

# # set_gate!(qc, 2, 1, "CNOT")
# set_gate!(qc, 2, 2, "X")
# H(qc, 1)
# CNOT(qc, 1, 2)
# CNOT(qc, 1, 2)
# CNOT(qc, 3, 2)
# CNOT(qc, 3, 4)
# R(qc, 1, 3.141, 0.0, 0.0)
# H(qc, 4)
# CNOT(qc, 4, 1)
# H(qc, 4)
# print_circuit(qc)
# A = compute_layer(qc, 2)
# print_matrix(A)
# A = compute_circ(qc)
# print_matrix(A)
# compute_sv(qc, state)
# print(state)
# print_state(state)
# print(A)
# println("⨁123"[1])
# println("⨁123"[1] == "⨁")
# println("⨁123"[1] == '⨁')
# qc = nothing # deleting QuantumCircuit freeing memory

# println(("Allo le monde"))

#start of testing 

#Task 1 GHZ state > |00 ..0> --> (  |00 ..0> +|11 ..1>  ) /sqrt(2) note |Qn , ... Q2,Q1>
function GHZ(num_qubits::Int)
  

    println("inilizing step")
    state = initialize_statevector(num_qubits, 0)
    qc = QuantumCircuit(num_qubits)
    # println(state)
    H(qc, 1)
    for i in 2:num_qubits  #first check point
        CNOT(qc, 1, i)
    end
    # mem = @allocated begin
    #     qc
    # end
    # println("Memory allocated for QuantumCircuit: $mem bytes")

    # ops::Vector{Matrix{Complex{Float64}}} = []
    # @time for i in 1:num_qubits #second check point
    #     push!(ops, compute_layer(qc, i))
    # end
    
    # @time for i in 1:num_qubits #third check point
    #     state .= ops[i] * state
    # # end
    # print("compute_step:")
    # V =compute_circ(qc)
    # print_matrix(V)
    @time compute_sv(qc, state)

    # print_circuit(qc) #fourth check point
    # save_circuit_to_csv(qc, "GHZ(25).csv")
    # print_state(state)
    # plot_probabilities(state)
    # qc = nothing
    print_circuit(qc)
    # empty!(ops)
    # println(state)
    # print_state(state)
    qc = nothing
    GC.gc()
end

# struct sv_tensor
#     state::Vector{ComplexF64}
#     qc::QuantumCircuit
#     ops::Vector{Matrix{Complex{Float64}}}
# end
    

# Profile.clear()
GHZ(9)

# Profile.print(maxdepth=10, mincount=10)
# GHZ(10)
# sleep(10)
# exit()