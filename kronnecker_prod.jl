using LinearAlgebra

abstract type AbstractMatrixNode{T} end

struct TerminalNode{T} <: AbstractMatrixNode{T}
    value::T
    dimensions::Tuple{Int,Int}
end

struct InternalNode{T} <: AbstractMatrixNode{T}
    children::Vector{AbstractMatrixNode{T}}
    dimensions::Tuple{Int,Int}
end

struct MatrixDD{T}
    root::AbstractMatrixNode{T}
end

function zero_node(T::Type, dims)
    TerminalNode(zero(T), dims)
end

function scalar_node(value::T, dims) where T
    TerminalNode(value, dims)
end

is_terminal(node::TerminalNode) = true
is_terminal(node::InternalNode) = false
is_zero(node::TerminalNode) = iszero(node.value)
is_zero(node::InternalNode) = false

function MatrixDD(matrix::AbstractMatrix{T}) where T
    MatrixDD(build_from_dense(matrix, 1:size(matrix,1), 1:size(matrix,2)))
end

function build_from_dense(matrix, rows, cols)
    dims = (length(rows), length(cols))
    if dims == (1,1)
        return TerminalNode(matrix[rows.start, cols.start], dims)
    end
    submatrix = view(matrix, rows, cols)
    if all(iszero, submatrix)
        return zero_node(eltype(submatrix), dims)
    end
    first_val = submatrix[1,1]
    if all(==(first_val), submatrix)
        return scalar_node(first_val, dims)
    end
    mr = rows.start + div(length(rows),2) - 1
    mc = cols.start + div(length(cols),2) - 1
    children = [
        build_from_dense(matrix, rows.start:mr, cols.start:mc),
        build_from_dense(matrix, rows.start:mr, (mc+1):cols.stop),
        build_from_dense(matrix, (mr+1):rows.stop, cols.start:mc),
        build_from_dense(matrix, (mr+1):rows.stop, (mc+1):cols.stop)
    ]
    InternalNode(children, dims)
end

function multiply(a::MatrixDD, b::MatrixDD)
    a.root.dimensions[2] == b.root.dimensions[1] || throw(DimensionMismatch())
    MatrixDD(multiply_nodes(a.root, b.root))
end

function multiply_nodes(a::TerminalNode, b::TerminalNode)
    (is_zero(a) || is_zero(b)) && return zero_node(a.value, (a.dimensions[1], b.dimensions[2]))
    scalar_node(a.value * b.value, (a.dimensions[1], b.dimensions[2]))
end

function multiply_nodes(a::AbstractMatrixNode, b::AbstractMatrixNode)
    (is_zero(a) || is_zero(b)) && return zero_node(a.value, (a.dimensions[1], b.dimensions[2]))
    if is_terminal(a) && is_terminal(b)
        return multiply_nodes(a, b)
    end
    a_children = is_terminal(a) ? expand_terminal(a) : a.children
    b_children = is_terminal(b) ? expand_terminal(b) : b.children
    result_dims = (a.dimensions[1], b.dimensions[2])
    result_children = Vector{AbstractMatrixNode}(undef, 4)
    for i in 1:2, j in 1:2
        sum_node = zero_node(eltype(a_children[1].value), (
            div(result_dims[1],2) + (i==2 ? rem(result_dims[1],2) : 0,
            div(result_dims[2],2) + (j==2 ? rem(result_dims[2],2) : 0))
        for k in 1:2
            prod = multiply_nodes(a_children[(i-1)*2 + k], b_children[(k-1)*2 + j])
            sum_node = add_nodes(sum_node, prod)
        end
        result_children[(i-1)*2 + j] = sum_node
    end
    collapse_node(InternalNode(result_children, result_dims))
end

function kronecker(a::MatrixDD, b::MatrixDD)
    MatrixDD(kronecker_nodes(a.root, b.root))
end

function kronecker_nodes(a::TerminalNode, b::TerminalNode)
    scalar_node(a.value * b.value, (a.dimensions[1]*b.dimensions[1], a.dimensions[2]*b.dimensions[2]))
end

function kronecker_nodes(a::AbstractMatrixNode, b::AbstractMatrixNode)
    if is_terminal(a) && is_terminal(b)
        return kronecker_nodes(a, b)
    end
    a_children = is_terminal(a) ? expand_terminal(a) : a.children
    b_children = is_terminal(b) ? expand_terminal(b) : b.children
    result_dims = (a.dimensions[1]*b.dimensions[1], a.dimensions[2]*b.dimensions[2])
    result_children = AbstractMatrixNode[]
    for a_child in a_children
        for b_child in b_children
            push!(result_children, kronecker_nodes(a_child, b_child))
        end
    end
    collapse_node(InternalNode(result_children, result_dims))
end

function multiply_vector(matrix::MatrixDD, vector)
    matrix.root.dimensions[2] == length(vector) || throw(DimensionMismatch())
    multiply_node_vector(matrix.root, vector)
end

function multiply_node_vector(node::TerminalNode, vector)
    fill(node.value * sum(vector), node.dimensions[1])
end

function multiply_node_vector(node::InternalNode, vector)
    mid = div(length(vector),2)
    v1 = @view vector[1:mid]
    v2 = @view vector[mid+1:end]
    top = multiply_node_vector(node.children[1], v1) .+ multiply_node_vector(node.children[2], v2)
    bottom = multiply_node_vector(node.children[3], v1) .+ multiply_node_vector(node.children[4], v2)
    vcat(top, bottom)
end

function expand_terminal(node::TerminalNode)
    c1 = scalar_node(node.value, (div(node.dimensions[1],2), div(node.dimensions[2],2)))
    c2 = scalar_node(node.value, (div(node.dimensions[1],2), node.dimensions[2] - div(node.dimensions[2],2)))
    c3 = scalar_node(node.value, (node.dimensions[1] - div(node.dimensions[1],2), div(node.dimensions[2],2)))
    c4 = scalar_node(node.value, (node.dimensions[1] - div(node.dimensions[1],2), node.dimensions[2] - div(node.dimensions[2],2)))
    [c1, c2, c3, c4]
end

function add_nodes(a::TerminalNode, b::TerminalNode)
    scalar_node(a.value + b.value, a.dimensions)
end

function add_nodes(a::AbstractMatrixNode, b::AbstractMatrixNode)
    if is_terminal(a) && is_terminal(b)
        return add_nodes(a, b)
    end
    a_children = is_terminal(a) ? expand_terminal(a) : a.children
    b_children = is_terminal(b) ? expand_terminal(b) : b.children
    result_children = [add_nodes(a_children[i], b_children[i]) for i in 1:4]
    collapse_node(InternalNode(result_children, a.dimensions))
end

function collapse_node(node::InternalNode)
    all_child_values = [c.value for c in node.children if is_terminal(c)]
    if length(all_child_values) == 4 && all(==(all_child_values[1]), all_child_values)
        return scalar_node(all_child_values[1], node.dimensions)
    end
    node
end

Base.:*(a::MatrixDD, b::MatrixDD) = multiply(a, b)
Base.:*(a::MatrixDD, v::AbstractVector) = multiply_vector(a, v)
Base.:kron(a::MatrixDD, b::MatrixDD) = kronecker(a, b)
