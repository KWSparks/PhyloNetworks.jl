using PhyloNetworks
import PhyloNetworks.getIndex
import PhyloNetworks.getOtherNode
import PhyloNetworks.addBL

net = readTopology("(((A:4.0,(B:1.0)#H1:1.1::0.9):0.5,(C:0.6,#H1:1.0::0.1):1.0):3.0,D:5.0);")
tips = Dict("A" => 0, "B" => 0, "C" => 1, "D" => 1)

"""
    tree_corelikelihood(tips, tree, mlik)

Calculate likelihood of discrete character states on a phylogenetic network given
starting states.

"""
function tree_corelikelihood(tree::HybridNetwork, tips::Dict{String,Int64}, mlik::Array{Float64}(length(net.node),size(k)))    
    for n in tree.node
        if n.leaf
            if !keys(n.name) #if tip name is not a key in dictionary
                mlik[:, n] = 0.0 #missing data: probability=1, log(1)=0.0
            end
        end
        if n.root
            lt = Array{length(k), 1}
            fill!(lt, log(1/k))
        else #get the parent edge index
            lt = logtrans[:,:,e]
        end
        if n.leaf
            if !n.root
                mlik[:, n] = lt[:,σ]
            else
                mlik[:, n] = 
            end
        end
        for i,j in mlik[i,j]
            mlik[i,n] = logsumexp(mlik[i,n], lt[i,j]+mlik[j,c])
            if i>=2 && n.root
                break
            elseif n.root
                return mlik[1,n]
            end
        end 
    end
end

"""
    network_corelikelihood(tips, mod, net)

Calculate likelihood for discrete characters on a network using 

# Examples

"""

network_corelikelihood(tips::Dict{String,Int64}, mod::TraitSubstitutionModel, trees::Array{HybridNetwork}, ltw)
    ll = 0
    f(t) = tree_corelikelihood(trees[t],view(mlik,:,:,t))
    pmap(f, 1:length(trees))
    for t in trees
        if t=1
            ll = mlik[1, root, t] + ltw[1]
        else
            ll = logsumexp(ll, mlik[1, root, t] + ltw[t])
    return ll
end

"""
    network_calculatelikelihood()

Calculate likelihood of discrete character states on a reticulate network given.

# Examples
"""

function network_calculatelikelihood(tips::Dict{Int64,Set{T}}, mod::TraitSubstitutionModel, net::HybridNetwork) # fixme: replace TraitSubstitutionModel with SM
    trees = displayedTrees(net, 0.0)
    k = nStates(mod)
    #both mlik and logtrans should be 3-d arrays
    #mlik[i,n,t] = log P{data below n in tree t given state i above n}
    #re-number edges to be consective, positive numbers; check edges are positive numbers
    #keep track of largest edge number
    #initialize 3d array: Array{Float64}((i,j,e))
    mlik = Array{Array{Float64, 1}(size(net.node))}(size(k))
    #logtrans[i,j,e]; i = start_state, j = end_state, e = edge.number
    logtrans = Array{Matrix{Float64}(size(k),size(k))}(length(net.edge))
    log_tree_weights = Array{Float64,1}
    ltw = 0.0
    #Step 1
    for tree in trees
        for e in tree.edge
            if e.gamma != 1.0
                ltw += log(e.gamma)
            end
        end
        push!(log_tree_weights,ltw)
        ltw = 0.0
    end
    #Step 2
    for edge in net.edge
        logtrans[:,:,edge.number] = log(P(mod,edge.length)) # element-wise
    end
    #Step 3
    network_likelihoodcore(tips,mod,net)
    #fixme: add optimization routine
    return
end