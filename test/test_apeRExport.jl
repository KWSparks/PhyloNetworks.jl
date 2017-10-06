R"library(ape)"

info("Testing `apeRExport` function")
#checking for absence of errors, not testing for correctness

#on a tree with some edge lengths missing
s = "(A,(B:1.0,(C:1.0,D:1.0):1.0):1.0);"
net = readTopology("(A,(B:1.0,(C:1.0,D:1.0):1.0):1.0);")
tree1 = apeRExport(net)
R"tree2 = read.evonet(text = $s)"
@test convert(Bool, R"isTRUE(all.equal(tree2,$tree1))")
        
#on a network with some missing gamma values
s = "(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:1.0):1.0):1.0,D:1.0);"
net = readTopology("(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:1.0):1.0):1.0,D:1.0);") 
tree1 = apeRExport(net)
R"tree2 = read.evonet(text = $s)"
@test convert(Bool, R"isTRUE(all.equal(tree2,$tree1))")

#on a network with mainTree=true; minor hybrid edge length missing
s = "(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:::0.1):1.0):1.0,D:1.0);"
net = readTopology("(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:::0.1):1.0):1.0,D:1.0);") 
tree1 = apeRExport(net; mainTree=true)
R"tree2 = read.evonet(text = $s)"
@test convert(Bool, R"isTRUE(all.equal(tree2,$tree1))")

#on a network with useEdgeLength=false 
net = readTopology("(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:1.0::0.1):1.0):1.0,D:1.0);")
tree1 = apeRExport(net; useEdgeLength=false)
R"tree2 = read.evonet(text = $s)"
@test convert(Bool, R"isTRUE(all.equal(tree2,$tree1))")

#on a networks with two hybridizations
net = readTopology("(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,(((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);")
s = "(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,(((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);"
tree1 = apeRExport(net)
R"tree2 = read.evonet(text = $s)"
@test convert(Bool, R"isTRUE(all.equal(tree2,$tree1))")

info("Testing `sexp` function")
#checking for absence of errors, not testing for correctness

#on a tree with some edge lengths missing
s = "(A,(B:1.0,(C:1.0,D:1.0):1.0):1.0);"
net = readTopology("(A,(B:1.0,(C:1.0,D:1.0):1.0):1.0);")
tree1 = sexp(net)
R"tree2 = read.evonet(text = $s)"
@test convert(Bool, R"isTRUE(all.equal(tree2,$tree1))")
        
#on a network with some missing gamma values
s = "(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:1.0):1.0):1.0,D:1.0);"
net = readTopology("(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:1.0):1.0):1.0,D:1.0);") 
tree1 = sexp(net)
R"tree2 = read.evonet(text = $s)"
@test convert(Bool, R"isTRUE(all.equal(tree2,$tree1))")

#on a network with mainTree=true; minor hybrid edge length missing
s = "(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:::0.1):1.0):1.0,D:1.0);"
net = readTopology("(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:::0.1):1.0):1.0,D:1.0);") 
tree1 = sexp(net; mainTree=true)
R"tree2 = read.evonet(text = $s)"
@test convert(Bool, R"isTRUE(all.equal(tree2,$tree1))")

#on a network with useEdgeLength=false 
net = readTopology("(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:1.0::0.1):1.0):1.0,D:1.0);")
tree1 = sexp(net; useEdgeLength=false)
R"tree2 = read.evonet(text = $s)"
@test convert(Bool, R"isTRUE(all.equal(tree2,$tree1))")

#on a networks with two hybridizations
net = readTopology("(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,(((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);")
s = "(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,(((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);"
tree1 = sexp(net)
R"tree2 = read.evonet(text = $s)"
@test convert(Bool, R"isTRUE(all.equal(tree2,$tree1))")
