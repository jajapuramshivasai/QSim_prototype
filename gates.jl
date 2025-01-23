using LinearAlgebra

function complex_array(n::Int, m::Int)::Vector{ComplexF64}
    size = 2^n
    result = zeros(ComplexF64, size)
    result[m+1] = 1.0
    return result
end




function multiply_and_update(arr::Vector{ComplexF64},  Unitary::Matrix{ComplexF64})
    # Check if the dimensions are compatible for multiplication
    if size(Unitary, 2) != length(arr)
        throw(ArgumentError("Matrix columns must match vector length for multiplication."))
    end
    
    # Perform the multiplication and update the vector
    arr .= Unitary * arr
end
function identity_matrix(qubits::Int)::Matrix{ComplexF64}
    n = 1<<qubits
    return Matrix{ComplexF64}(I, n, n)
end








function Hadamard(arr::Vector{ComplexF64}, target::Int)
    H::Matrix{ComplexF64} =  1/sqrt(2) * [1 1; 1 -1]
    qubits::Int = log2(length(arr))
    # arr .= kron(identity_matrix(qubit-1), kron(H, identity_matrix(n-qubit))) * arr
    l::Int = target - 1 
    r::Int= qubits - target

    # println( kron(identity_matrix(l),kron(H,identity_matrix(r))))
    arr .= kron(identity_matrix(r),kron(H,identity_matrix(l))) * arr
end

function Pauli_X(arr::Vector{ComplexF64}, target::Int)
    X::Matrix{ComplexF64} = [0 1; 1 0]
    qubits::Int = log2(length(arr))
    # arr .= kron(identity_matrix(qubit-1), kron(H, identity_matrix(n-qubit))) * arr
    l::Int = target - 1 
    r::Int= qubits - target

    # println( kron(identity_matrix(l),kron(H,identity_matrix(r))))
    arr .= kron(identity_matrix(r),kron(X,identity_matrix(l))) * arr
end

function Pauli_Y(arr::Vector{ComplexF64}, target::Int)
    Y::Matrix{ComplexF64} = [0 -im; im 0]
    qubits::Int = log2(length(arr))
    # arr .= kron(identity_matrix(qubit-1), kron(H, identity_matrix(n-qubit))) * arr
    l::Int = target - 1 
    r::Int= qubits - target

    # println( kron(identity_matrix(l),kron(H,identity_matrix(r))))
    arr .= kron(identity_matrix(r),kron(Y,identity_matrix(l))) * arr
end

function Pauli_Z(arr::Vector{ComplexF64}, target::Int)
    Z::Matrix{ComplexF64} = [1 0; 0 -1]
    qubits::Int = log2(length(arr))
    # arr .= kron(identity_matrix(qubit-1), kron(H, identity_matrix(n-qubit))) * arr
    l::Int = target - 1 
    r::Int= qubits - target

    # println( kron(identity_matrix(l),kron(H,identity_matrix(r))))
    arr .= kron(identity_matrix(r),kron(Z,identity_matrix(l))) * arr
end


function Unitary(arr::Vector{ComplexF64}, target::Int, U::Matrix{ComplexF64})
    qubits::Int = log2(length(arr))
    r::Int = target - 1 
    l::Int= qubits - target
    arr .= kron(identity_matrix(l), kron(U, identity_matrix(r))) * arr
    
end

function CNOT(arr::Vector{ComplexF64}, control::Int, target::Int)
    qubits::Int = log2(length(arr))
    p0::Matrix{ComplexF64} = [1 0; 0 0]
    p1::Matrix{ComplexF64} = [0 0; 0 1]
    X::Matrix{ComplexF64} = [0 1; 1 0]
    I::Matrix{ComplexF64} = identity_matrix(1)

    control, target = target, control

    if control < target
        l = control - 1 
        m = target - control - 1
        r = qubits - target

        # println("right: ", r)
        # println("middle: ", m)  
        # println("left: ", l)
        # println([r,"C",m,"T",l])


    # U1 = kron((identity_matrix(l), p0) ,kron(identity_matrix(m) ,kron(I, identity_matrix(r))))
    U1 = kron(kron(identity_matrix(l), p0), 
          kron(identity_matrix(m), 
               kron(I, identity_matrix(r))))

    # U2 = kron((identity_matrix(l), p1) ,kron(identity_matrix(m) ,kron(X, identity_matrix(r))))
    U2 = kron(kron(identity_matrix(l), p1), 
          kron(identity_matrix(m), 
               kron(X, identity_matrix(r))))

    U = U1 + U2
    println(U)
    arr .= U * arr
    else
        control, target = target, control
        l = control - 1 
        m = target - control - 1
        r = qubits - target

        # println("right: ", r)
        # println("middle: ", m)  
        # println("left: ", l)
        # println([r,"T",m,"C",l])

    # U1 = kron((identity_matrix(r), p0) ,kron(identity_matrix(m) ,kron(I, identity_matrix(l))))
    U1 = kron(kron(identity_matrix(l), I), 
    kron(identity_matrix(m), 
         kron(p0, identity_matrix(r))))
    # U2 = kron((identity_matrix(r), p1) ,kron(identity_matrix(m) ,kron(X, identity_matrix(l))))
    U2 = kron(kron(identity_matrix(l), X), 
    kron(identity_matrix(m), 
         kron(p1, identity_matrix(r))))
    U = U1 + U2
    # println(U)
    arr .= U * arr
    end


    # main logic
    # l::Int = control - 1 
    # m::Int = target - control - 1
    # r::Int= qubits - target

    #debugging
    # println("right: ", r)
    # println("middle: ", m)  
    # println("left: ", l)
    

    # U1 = kron((identity_matrix(l), p0) ,kron(identity_matrix(m) ,kron(I, identity_matrix(r))))
    # U2 = kron((identity_matrix(l), p1) ,kron(identity_matrix(m) ,kron(X, identity_matrix(r))))
    # println(U1)
end





# Experimental code for Circuit distribution and evaluation

struct QuantumCircuit
    gates::Vector{Vector{String}}
    num_qubits::Int
    
    function QuantumCircuit(num_qubits::Int)
        new(Vector{Vector{String}}(), num_qubits)
    end
end

function add_layer!(circuit::QuantumCircuit)
    push!(circuit.gates, fill("I", circuit.num_qubits))
end

function set_gate!(circuit::QuantumCircuit, layer::Int, qubit::Int, gate::String)
    if layer > length(circuit.gates)
        for _ in 1:(layer - length(circuit.gates))
            add_layer!(circuit)
        end
    end
    circuit.gates[layer][qubit] = gate
end

function print_circuit(circuit::QuantumCircuit)
    for qubit in 1:circuit.num_qubits
        print("Q$qubit: ")
        for layer in circuit.gates
            print(layer[qubit], " ")
        end
        println()
    end
end
# end of experimental code

function measure_probabilities(state_vector::Vector{ComplexF64})::Vector{Float64}
    # Calculate probabilities by taking the absolute square of amplitudes
    return abs2.(state_vector)
end

# using Plots

# function plot_probabilities(state_vector::Vector{ComplexF64})
#     probs = measure_probabilities(state_vector)
#     n = length(probs)
#     bar(0:(n-1), probs, 
#         title="State Probabilities",
#         xlabel="State",
#         ylabel="Probability",
#         legend=false)
# end

function print_state_representation(state_vector::Vector{ComplexF64})
    n = Int(log2(length(state_vector)))
    println("State representation:")
    for i in 0:(length(state_vector)-1)
        amplitude = state_vector[i+1]
        if abs(amplitude) > 1e-10  # Only print non-zero amplitudes
            binary = lpad(string(i, base=2), n, '0')
            println("|$binary⟩: $amplitude")
        end
    end
end





# Main
println("Allo le Monde!")
for i in 1:10
    print("++++++++++++++")
end


"""
Notation ++> md"|0(4) 0(3) 0(2) 0(1)>"

for example

    lets initialize two qubits in       |00>

    after applying Pauli_X on qubit 1
    we get                              |01>

    after applying CNOT with control 
    on qubit 1 and target on qubit 2 
    we get                              |11>

"""


println("initializing vector")
v = complex_array(2,0)
# println(v)
print_state_representation(v)

println("applying Pauli_X")
# Pauli_X(v, 2)
# print_state_representation(v)
Pauli_X(v, 1)
print_state_representation(v)
# println(v)

# println("applying matrix")
# multiply_and_update(v, Pauli_X)
# println(v)

# println("initializing Identity matrix")
# I2 = identity_matrix(2)
# println(I2)

# println("applying Hadamard")
# Hadamard(v, 1)
# println(v)

# println("applying Unitary")
# Pauli_X::Matrix{ComplexF64} = [0 1; 1 0]
# Unitary(v, 1, Pauli_X)
# println(v)

println("applying CNOT")
CNOT(v, 1, 2)
# println(v)
print_state_representation(v)

# measuring and plotting the probabilities
# plot_probabilities(v)

#exit routine
println("Task executed successfully killing process in 5 seconds")
sleep(10)
exit()
