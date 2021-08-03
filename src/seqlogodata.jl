### Coordinates for ACGT
## need to replace with something betterm, uses matplotlib textpath
function motif_letter_data(; s=10, weight="bold")
       
       fp = PyPlot.matplotlib.font_manager.FontProperties(family="arial", weight=weight)

       a = PyPlot.matplotlib.textpath.TextPath((0, 0), "A", size=1, prop=fp)
       c = PyPlot.matplotlib.textpath.TextPath((0, 0), "C", size=s, prop=fp) ## setting size as 400 gives > 100 interp points for bezier curve
       g = PyPlot.matplotlib.textpath.TextPath((0, 0), "G", size=s, prop=fp)
       t = PyPlot.matplotlib.textpath.TextPath((0, 0), "T", size=1, prop=fp)
              
       adata =  lettertuple(a.to_polygons())
       cdata  = lettertuple([PyPlot.matplotlib.patches.PathPatch(c).get_verts()])
       gdata  = lettertuple([PyPlot.matplotlib.patches.PathPatch(g).get_verts()])
       tdata  = lettertuple(t.to_polygons())

       _agct = Dict(DNA_A => adata, DNA_C => cdata, DNA_G => gdata, DNA_T => tdata)
       agct = Dict(b => normletter(l) for (b, l) in _agct)

       agct
end
