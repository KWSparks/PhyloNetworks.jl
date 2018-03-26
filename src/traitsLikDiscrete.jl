using PhyloNetworks
import PhyloNetworks.getIndex
import PhyloNetworks.getOtherNode
import PhyloNetworks.addBL

net = readTopology("(((A:4.0,(B:1.0)#H1:1.1::0.9):0.5,(C:0.6,#H1:1.0::0.1):1.0):3.0,D:5.0);")
tips = Dict("A" => 0, "B" => 0, "C" => 1, "D" => 1)

"""
    discrete_tree_corelikelihood(tree, tips, logtrans)

Calculate likelihood of discrete character states on a phylogenetic network given
starting states.

"""
function discrete_tree_corelikelihood(tree::HybridNetwork, tips::Dict{String,Int64},
    logtrans::AbstractArray, fixit
    forwardlik = Array{Array{Float64, 1}(size(net.node))}(size(k))
    directlik = Array{Array{Float64, 1}(size(net.edge))}(size(k))
    backwardlik = Array{Array{Float64, 1}(size(net.node))}(size(k)))
    for i in reverse(1:length(tree.nodes_changed)) # post-order
        n = tree.nodes_changed[i]
        if n.leaf
            if !keys(n.name) #if tip name is not a key in dictionary
                for i in 1:k
                    forwardlik[i,n] = 0.0 
                end
            else
                for i in 1:k
                    if i = tips(n.name)
                        forwardlik[i,n]
                    else
                        forwardlik[i,n] = -Inf64
                    end
                end
            end
        elseif n.root
            for i in 1:k
                tmp = logprior[i] + forwardlik[i,n]
                if i=1
                    loglik = tmp
                elseif i > 1
                    loglik = logsumexp(loglik, tmp)
                else
                    for e in n.edge
                        pe = e.isChild1 ? 1 : 2
                    end
                    lt = logtrans[:,:,e]
                    directlik[:,e] = logtrans[:,1,e] + forwardlik[1,n]
                    for i in 1:k
                        for j in 2:k
                            tmp = logtrans[i,j,e] + forwardlik[j,n]
                            directlik[i,e] = logsumexp(directlik[i,e],tmp)
                        end
                    end
                end
            end
        else
            children = []
            for e in n.edge
                if e.node[e.isChild1 ? 1 : 2] == node continue; end
                # excluded parent edges only: assuming tree here
                push!(children, e)
            end
            e1 = children[1]
            e2 = children[2]
            forwardlik[:,n] = directlik[:,e1] + directlik[:,e2]
        ende.node[e.isChild1 ? 1 : 2]
    end
    for n in net.nodes_changed
        if n.root
            logprior = log(backwardlik[:,n])
        else
            pn = n.isChild1 ? 1 : 2
            for e in n.edge
                pe = e.isChild1 ? 1 : 2
            end 
            for i in 1:k
                for j in 1:k
                    tmp = backwardlik[j,pn] +logtrans[j,i,pn] + # Fixit: sum child of pn: directlik[j,e]
                    if j=1
                        backwardlik[i,n] = tmp
                    elseif j > 1
                        backwardlik[i,n] = logsumexp(backwardlik[i,n],tmp)
                    end
                end
            end
        end
    end
    return loglik
    #Fixit: ancestral state reconstruction at node n
end

"""
    discrete_corelikelihood(tips, mod, trees, ltw,
        logtrans,forwardlik,directlik,backwardlik)

Calculate likelihood for discrete characters on a network,
using fixed model parameters

# Examples

"""

function discrete_corelikelihood(tips::Dict{String,Int64}, mod::TraitSubstitutionModel,
    trees::Array{HybridNetwork}, ltw::AbstractVector, logtrans::AbstractArray,
    forwardlik::AbstractArray, directlik::AbstractArray, backwardlik::AbstractArray)
    f(t) = discrete_tree_corelikelihood(trees[t],logtrans,view(forwardlik,:,:,t),
             view(directlik,:,:,t),view(backwardlik,:,:,t))
    ll = pmap(f, 1:length(trees)) # ll = loglikelihood given each tree
    @show ll
    res = ll[1] + ltw[1] # result: loglikelihood given the network
    @show res
    for t in 2:length(trees)
        res = logsumexp(res, ll[t] + ltw[t])
    end
    return res
end

"""
    discrete_optimlikelihood(tips, mod, net)

Calculate likelihood of discrete character states on a reticulate network.

"""

function discrete_optimlikelihood(tips::Dict{Int64,Set{T}}, mod::SM, net::HybridNetwork)
    trees = displayedTrees(net, 0.0)
    ntrees = length(trees)
    k = nStates(mod)
    #both mlik and logtrans should be 3-d arrays
    #mlik[i,n,t] = log P{data below n in tree t given state i above n}
    # fixit: re-number edges to be consective, positive numbers; check edges are positive numbers
    #keep track of largest edge number
    #initialize 3d array: Array{Float64}((i,j,e))
    # fixit: later change the arrays below into SharedArray
    forwardlik = Array{Float64}(k,size(net.node),ntrees)
    directlik  = Array{Float64}(k,size(net.edge),ntrees)
    backwardlik= Array{Float64}(k,size(net.node),ntrees)
    #logtrans[i,j,e]; i = start_state, j = end_state, e = edge.number
    log_tree_weights = Array{Float64,1}
    #Step 1
    for tree in trees
        ltw = 0.0
        for e in tree.edge
            if e.gamma != 1.0
                ltw += log(e.gamma)
            end
        end
        push!(log_tree_weights,ltw)
    end
    #Step 2
    logtrans = Array{Float64}(k,k,length(net.edge))
    for edge in net.edge
        logtrans[:,:,edge.number] = log(P(mod,edge.length)) # element-wise
    end
    #Step 3
    discrete_corelikelihood(tips,mod,trees,ltw,logtrans,forwardlik,directlik,backwardlik)
    #fixme: add optimization routine
    return
end


