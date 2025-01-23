

"""
To Do:


>implement intermediate representation of circuit
>pack the quantum gates as far as possible
>>Example : <{I,X},{H,I}> -> <{H,I}> #basically eliminate as many I as possible
>>Example : <{I,I,,I,C,-,X},{C,-,X,I,I,I}> -> <{C,-,X},{C,-,X}> #instead of comma we can use empty string or space

>convert string matrix into array of stack of complex matrices representing quantum gates 
>>Example : {I,I,C0,X1} -> {I(2),CNOT(2,1)} -> {[[1 0 0 0],[0 1 0 0],[0 0 1 0], [0 0 0 1]], [[1 0 0 0],[0 0 0 1],[0 0 1 0], [0 1 0 0]]}
>>Example : {CNOT(3,1)} #complete this line

>convert array of stack of complex matrices into array of refrence to the unitaries of each circuit layer efficiently represented in desicion diagram data structure
>>Example : #complete this line 

>Apply the circuit layer unitaries onto state vector iteratively

!old module : "gates.jl"
>now with above structure we can implement tensor networks and quantum circuits efficiently
>add structure QuantumStateVector initiall storing state vector as array of individual qubit state_vectors
>>ex v = QuantumStateVector(qubits = 2);print(v) => {[1,0],[1,0]}
>when two qubit gate is applied contrach the graph and update the state vector
>>ex Pauli_X(v, 1);print(v) => {[0,1],[1,0]} 
>>ex CNOT(v, 1, 2);print(v) => {[0,0,0,1]}
>>ex u = {[0,1],[0,1],[0,1],[0.71,0.71]};CNOT(u,1,2);print(u) => {[0,1],[0,1],[0.71,0,0,0.71]}
>use tensor networks to compute sub circuits as single unitary

!new module : "gradient"

>implement multi circuit jobs that take multiple circuits and applies one after another or combine and compute the unitary of entire sub-circuits
>implement parameterized sub circuits
>implement gradient descent on the parameters of the sub-circuits using parameter shift rule
>implement gradient descent on the parameters of the sub-circuits using finite difference method
>implement gradient descent on the parameters of the sub-circuits using backpropagation

>intigrate with the flux.jl module to implement the quantum neural network

"""
struct Circuit
    data::Array{String, 2}
    stacks::Vector{Vector{Matrix{Complex{Float64}}}}
end

function Circuit(rows::Int, cols::Int)
    return Circuit(fill("", rows, cols), Vector{Vector{Matrix{Complex{Float64}}}}())
end

function add_stack!(circuit::Circuit, stack::Vector{Matrix{Complex{Float64}}})
    push!(circuit.stacks, stack)
end

function get_stack(circuit::Circuit, index::Int)
    return circuit.stacks[index]
end
function iterate_stack(circuit::Circuit, index::Int)
    stack = get_stack(circuit, index)
    for matrix in stack
        println(matrix)
    end
end

# Example to test Circuit object
function test_circuit()
    # Create a Circuit object with 2 rows and 3 columns
    circuit = Circuit(2, 3)
    
    # Create a stack of matrices
    matrix1 = [Complex{Float64}(1.0, 2.0)  Complex{Float64}(3.0, 4.0); Complex{Float64}(5.0, 6.0) Complex{Float64}(7.0, 8.0)]
    matrix2 = [Complex{Float64}(9.0, 10.0) Complex{Float64}(11.0, 12.0); Complex{Float64}(13.0, 14.0) Complex{Float64}(15.0, 16.0)]
    stack = [matrix1, matrix2]
    CNOT::Matrix{Complex{Float64}} = [0 0 0 1; 0 1 0 0; 0 0 1 0; 1 0 0 0]
    # Add the stack to the circuit
    add_stack!(circuit, stack)
    add_stack!(circuit, [CNOT])
    
    # Retrieve and print the stack from the circuit
    retrieved_stack = get_stack(circuit, 1)
    println("Retrieved stack:")
    for matrix in retrieved_stack
        println(matrix)
    end
    
    # Iterate and print each matrix in the stack
    println("Iterating stack:")
    iterate_stack(circuit, 2)
end

# Run the test
test_circuit()

sleep(5)
exit()
