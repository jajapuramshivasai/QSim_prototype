using LinearAlgebra
using SparseArrays

abstract type AbstractDDNode{T} end

struct TerminalNode{T} <: AbstractDDNode{T}
    value::Union{T, Nothing}
    dimensions::Tuple{Int,Int}
    is_zero::Bool
end

struct InternalNode{T} <: AbstractDDNode{T}
    children::Vector{AbstractDDNode{T}}
    split_axis::Symbol  # :row, :col, :both
    dimensions::Tuple{Int,Int}
    hash::UInt
end

struct UnitaryDD{T} <: AbstractMatrix{T}
    root::AbstractDDNode{T}
    sparse_map::Dict{UInt,AbstractDDNode{T}}
end

function Base.size(u::UnitaryDD)
    u.root.dimensions
end

function Base.getindex(u::UnitaryDD, i::Int, j::Int)
    function _getindex(node, i, j)
        if node isa TerminalNode
            return node.is_zero ? zero(eltype(u)) : node.value
        end
        
        r, c = node.dimensions
        if node.split_axis == :row
            child_idx = i <= div(r,2) ? 1 : 2
            new_i = child_idx == 1 ? i : i - div(r,2)
            _getindex(node.children[child_idx], new_i, j)
        elseif node.split_axis == :col
            child_idx = j <= div(c,2) ? 1 : 2
            new_j = child_idx == 1 ? j : j - div(c,2)
            _getindex(node.children[child_idx], i, new_j)
        else
            row_child = i <= div(r,2) ? 1 : 3
            col_child = j <= div(c,2) ? 1 : 2
            child_idx = row_child + (col_child - 1)
            new_i = child_idx in [1,2] ? i : i - div(r,2)
            new_j = child_idx in [1,3] ? j : j - div(c,2)
            _getindex(node.children[child_idx], new_i, new_j)
        end
    end
    _getindex(u.root, i, j)
end

function UnitaryDD(matrix::AbstractMatrix{T}; ϵ=1e-9) where T
    root = build_node(matrix, 1:size(matrix,1), 1:size(matrix,2), ϵ)
    UnitaryDD(root, Dict{UInt,AbstractDDNode{T}}())
end

function build_node(mat, rows, cols, ϵ)
    r_len, c_len = length(rows), length(cols)
    sub = view(mat, rows, cols)
    
    is_unitary = norm(sub'*sub - I) < ϵ
    all_zero = all(x->abs(x)<ϵ, sub)
    const_val = all(x≈sub[1] for x in sub) ? sub[1] : nothing
    
    if const_val !== nothing
        return TerminalNode(const_val, (r_len,c_len), false)
    elseif all_zero
        return TerminalNode(nothing, (r_len,c_len), true)
    end
    
    split_row = r_len > 1 && (c_len == 1 || r_len >= c_len)
    split_col = c_len > 1 && !split_row
    
    children = AbstractDDNode{T}[]
    if split_row
        mid = div(r_len,2)
        push!(children, build_node(mat, rows[1:mid], cols, ϵ))
        push!(children, build_node(mat, rows[mid+1:end], cols, ϵ))
        axis = :row
    elseif split_col
        mid = div(c_len,2)
        push!(children, build_node(mat, rows, cols[1:mid], ϵ))
        push!(children, build_node(mat, rows, cols[mid+1:end], ϵ))
        axis = :col
    else
        r_mid = div(r_len,2)
        c_mid = div(c_len,2)
        push!(children, build_node(mat, rows[1:r_mid], cols[1:c_mid], ϵ))
        push!(children, build_node(mat, rows[1:r_mid], cols[c_mid+1:end], ϵ))
        push!(children, build_node(mat, rows[r_mid+1:end], cols[1:c_mid], ϵ))
        push!(children, build_node(mat, rows[r_mid+1:end], cols[c_mid+1:end], ϵ))
        axis = :both
    end
    
    node = InternalNode(children, axis, (r_len,c_len), hash(children))
    collapse!(node)
end

function collapse!(node::InternalNode{T}) where T
    unique_children = unique(node.children)
    if length(unique_children) == 1
        child = unique_children[1]
        if child isa TerminalNode
            return child
        end
    end
    
    all_zero = true
    for c in node.children
        if !(c isa TerminalNode) || !c.is_zero
            all_zero = false
            break
        end
    end
    all_zero && return TerminalNode(nothing, node.dimensions, true)
    
    node
end

function Base.:*(a::UnitaryDD{T}, b::UnitaryDD{T}) where T
    size(a,2) == size(b,1) || throw(DimensionMismatch())
    UnitaryDD(_mul_root(a.root, b.root))
end

function _mul_root(a::AbstractDDNode{T}, b::AbstractDDNode{T}) where T
    a_split = a isa InternalNode ? a.split_axis : :none
    b_split = b isa InternalNode ? b.split_axis : :none
    
    if a isa TerminalNode && b isa TerminalNode
        val = a.is_zero || b.is_zero ? zero(T) : a.value * b.value
        return TerminalNode(val, (a.dimensions[1], b.dimensions[2]), a.is_zero | b.is_zero)
    end
    
    if a_split == :col && b_split == :row
        children = AbstractDDNode{T}[]
        for a_child in a.children
            for b_child in b.children
                push!(children, _mul_root(a_child, b_child))
            end
        end
        new_node = InternalNode(children, :both, (a.dimensions[1], b.dimensions[2]), hash(children))
        return collapse!(new_node)
    end
    
    error("Complex splitting pattern not implemented")
end

function tensor_product(a::UnitaryDD{T}, b::UnitaryDD{T}) where T
    UnitaryDD(_tensor_root(a.root, b.root))
end

function _tensor_root(a::AbstractDDNode{T}, b::AbstractDDNode{T}) where T
    if a isa TerminalNode
        a.is_zero && return TerminalNode(nothing, (a.dimensions[1]*b.dimensions[1], a.dimensions[2]*b.dimensions[2]), true)
        return _scale_node(b, a.value)
    end
    
    children = AbstractDDNode{T}[]
    for child in a.children
        push!(children, _tensor_root(child, b))
    end
    new_dims = (a.dimensions[1]*b.dimensions[1], a.dimensions[2]*b.dimensions[2])
    InternalNode(children, a.split_axis, new_dims, hash(children)) |> collapse!
end

function _scale_node(node::AbstractDDNode{T}, factor::T) where T
    node isa TerminalNode && return TerminalNode(node.value * factor, node.dimensions, false)
    InternalNode([_scale_node(c, factor) for c in node.children], node.split_axis, node.dimensions, hash(node.children)) |> collapse!
end

function print_dd(u::UnitaryDD; indent=0)
    function _print(node, level)
        prefix = "  "^level
        if node isa TerminalNode
            println(prefix, node.is_zero ? "𝟎" : "𝐔(≈$(round(node.value; digits=2))", 
                   " [$(node.dimensions[1])×$(node.dimensions[2])])")
        else
            println(prefix, "■ split:$(node.split_axis) [$(node.dimensions[1])×$(node.dimensions[2])]")
            for c in node.children
                _print(c, level+1)
            end
        end
    end
    _print(u.root, indent)
end

function Base.show(io::IO, u::UnitaryDD)
    print(io, "UnitaryDD$(size(u)) with ")
    nnz = count(!iszero, Matrix(u))
    print(io, "$nnz non-zero elements")
end

function contract!(u::UnitaryDD)
    function _contract(node)
        node isa TerminalNode && return node
        unique_children = unique(node.children)
        if length(unique_children) == 1
            return unique_children[1]
        end
        InternalNode([_contract(c) for c in node.children], node.split_axis, node.dimensions, hash(node.children))
    end
    u.root = _contract(u.root) |> collapse!
    u
end


# Example:
# id = UnitaryDD(Matrix(1.0I, 4,4))

# U = randn(4,4) + im*randn(4,4)
# U /= sqrt(tr(U'*U))
# ud = UnitaryDD(U)

# result = ud * id  # Matrix multiplication
# kron_prod = tensor_product(ud, id)  # Tensor product
# contract!(kron_prod)  # Optimize memory

# print_dd(result)
