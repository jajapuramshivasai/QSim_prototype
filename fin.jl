#issure:replace "I" with "-" #status:fixed
#issure:replace "⨁" with "⊕" or "⨁"

#implement HQNN and quantum transformer
using BenchmarkTools
using LinearAlgebra
using Plots

struct QuantumCircuit
    gates::Vector{Vector{String}}
    num_qubits::Int
    ops::Vector{Vector{Matrix{Complex{Float64}}}}
    function QuantumCircuit(num_qubits::Int)
        circuit = new(Vector{Vector{String}}(), num_qubits, Vector{Vector{Matrix{Complex{Float64}}}}(undef, num_qubits))
        add_layer!(circuit)
        return circuit
    end
end

function set_op!(circuit::QuantumCircuit, qubit::Int, op::Matrix{Complex{Float64}})
    if circuit.ops[qubit] === nothing
        circuit.ops[qubit] = Stack{Matrix{Complex{Float64}}}()
    end
    push!(circuit.ops[qubit], op)
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
            circuit.gates[layer+1][target] = "⨁$d"
            circuit.gates[layer+1][control] = "●$d"
            for i in control+1:target-1
                circuit.gates[layer+1][i] = "|"
            end
        else
            circuit.gates[layer][target] = "⨁$d"
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
            circuit.gates[layer+1][target] = "⨁$d"
            for i in target+1:control-1
                circuit.gates[layer+1][i] = "|"
            end
        else
            circuit.gates[layer][control] = "●$d"
            circuit.gates[layer][target] = "⨁$d"
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

function print_circuit(circuit::QuantumCircuit)
    for qubit in 1:circuit.num_qubits
        print("Q$qubit: ")
        for layer in circuit.gates
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
    for i in 0:(length(state_vector)-1)
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
    for i in 0:(length(state_vector)-1)
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
    while qubit <= circuit.num_qubits
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
    for layer in 1:length(circuit.gates)
        U = compute_layer(circuit, layer)
        state .= U * state
    end
    
end

function compute_circ(circuit::QuantumCircuit)
    A::Matrix{ComplexF64} = I(1<<circuit.num_qubits)
    for layer in 1:length(circuit.gates)
        U = compute_layer(circuit, layer)
        A .= U * A
    end
    return A
    
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
    state = initialize_statevector(num_qubits, 0)
    qc = QuantumCircuit(num_qubits)


    H(qc, 1)
    for i in 2:num_qubits
        CNOT(qc, 1, i)
    end
    print_circuit(qc)
    compute_sv(qc, state)
    print_state(state)
    # plot_probabilities(state)
    
end

@time GHZ(14)
