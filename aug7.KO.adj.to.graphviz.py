"""
digraph {
  graph [fontname = "helvetica", nodesep=0.300000, splines=ortho, ranksep=0.400000, rank=same];
  node [fontname = "helvetica", color="#000000", style=filled, fillcolor="#ffffff"];
  edge [fontname = "helvetica", penwidth=1];
  "154" -> "107"[color="#222222", penwidth=4.512063];
  "107" -> "343"[color="#4197c7", constraint=false, dir=none, penwidth=6.054689,style=dashed];
}
"""
import matrix_io as mio
import numpy as np

D = mio.load("KO.adj.matrix.tab", dtype=np.int)
M = D['M']
assert M.shape[0]==M.shape[1], M.shape
n = M.shape[0]
for i in xrange(n):
  for j in xrange(n):
    options = {}
    if M[i,j] == 1: # weak activator
      options.update({'color':'#339900'})
    elif M[i,j] == 2: # strong activator
      options.update({'color':'#00ff00'})
    elif M[i,j] == -1: # weak repressor
      options.update({'color':'orange', 'arrowhead':'tee'})
    elif M[i,j] == -2: # strong repressor
      options.update({'color':'red', 'arrowhead':'tee'})
    else:
      continue
    opts = ", ".join(['%s="%s"'%(k,v) for k,v in options.items()])
    print '"%s" -> "%s"[%s]' % (D['col_ids'][j], D['row_ids'][i], opts)
  
