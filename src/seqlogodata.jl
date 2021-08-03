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

function motif_letter_data_orig()
       a1 = [ (0.      , 0.      ),
       (3.345   , 8.589375),
       (5.180625, 8.589375),
              (8.619375, 0.      ),
              (6.733125, 0.      ),
              (5.983125, 1.951875),
              (2.548125, 1.951875),
              (1.839375, 0.      ),
              (0.      , 0.      )]

       a2  = [(5.42625 , 3.399375),
              (4.243125, 6.586875),
              (3.0825  , 3.399375)]

       a = [a1,  a2]

       c =   [( 6.369375,  3.1575  ),
              ( 8.05125 ,  2.625   ),
              ( 7.665   ,  1.21875 ),
              ( 6.765   ,  0.53625 ),
              ( 5.865   , -0.14625 ),
              ( 4.483125, -0.14625 ),
              ( 2.77125 , -0.14625 ),
              ( 1.66875 ,  1.021875),
              ( 0.568125,  2.191875),
              ( 0.568125,  4.21875 ),
              ( 0.568125,  6.36375 ),
              ( 1.674375,  7.54875 ),
              ( 2.7825  ,  8.735625),
              ( 4.588125,  8.735625),
              ( 6.165   ,  8.735625),
              ( 7.149375,  7.805625),
              ( 7.734375,  7.254375),
              ( 8.026875,  6.223125),
              ( 6.31125 ,  5.8125  ),
              ( 6.1575  ,  6.48    ),
              ( 5.67375 ,  6.86625 ),
              ( 5.191875,  7.254375),
              ( 4.5     ,  7.254375),
              ( 3.545625,  7.254375),
              ( 2.949375,  6.568125),
              ( 2.355   ,  5.88375 ),
              ( 2.355   ,  4.348125),
              ( 2.355   ,  2.71875 ),
              ( 2.94    ,  2.026875),
              ( 3.526875,  1.336875),
              ( 4.464375,  1.336875),
              ( 5.15625 ,  1.336875),
              ( 5.653125,  1.775625),
              ( 6.151875,  2.214375),
              ( 6.369375,  3.1575  )]

       c =   [( 0.53078125,  0.263125  ),
              ( 0.6709375 ,  0.21875   ),
              ( 0.63875   ,  0.1015625 ),
              ( 0.56375   ,  0.0446875 ),
              ( 0.48875   , -0.0121875 ),
              ( 0.37359375, -0.0121875 ),
              ( 0.2309375 , -0.0121875 ),
              ( 0.1390625 ,  0.08515625),
              ( 0.04734375,  0.18265625),
              ( 0.04734375,  0.3515625 ),
              ( 0.04734375,  0.5303125 ),
              ( 0.13953125,  0.6290625 ),
              ( 0.231875  ,  0.72796875),
              ( 0.38234375,  0.72796875),
              ( 0.51375   ,  0.72796875),
              ( 0.59578125,  0.65046875),
              ( 0.64453125,  0.60453125),
              ( 0.66890625,  0.51859375),
              ( 0.5259375 ,  0.484375  ),
              ( 0.513125  ,  0.54      ),
              ( 0.4728125 ,  0.5721875 ),
              ( 0.43265625,  0.60453125),
              ( 0.375     ,  0.60453125),
              ( 0.29546875,  0.60453125),
              ( 0.24578125,  0.54734375),
              ( 0.19625   ,  0.4903125 ),
              ( 0.19625   ,  0.36234375),
              ( 0.19625   ,  0.2265625 ),
              ( 0.245     ,  0.16890625),
              ( 0.29390625,  0.11140625),
              ( 0.37203125,  0.11140625),
              ( 0.4296875 ,  0.11140625),
              ( 0.47109375,  0.14796875),
              ( 0.51265625,  0.18453125),
              ( 0.53078125,  0.263125  )]

       g =   [( 0.40578125,  0.263125  ),
              ( 0.40578125,  0.38375   ),
              ( 0.71734375,  0.38375   ),
              ( 0.71734375,  0.09859375),
              ( 0.671875  ,  0.0546875 ),
              ( 0.585625  ,  0.02125   ),
              ( 0.49953125, -0.0121875 ),
              ( 0.41109375, -0.0121875 ),
              ( 0.29890625, -0.0121875 ),
              ( 0.2153125 ,  0.03484375),
              ( 0.131875  ,  0.08203125),
              ( 0.08984375,  0.1696875 ),
              ( 0.0478125 ,  0.25734375),
              ( 0.0478125 ,  0.3603125 ),
              ( 0.0478125 ,  0.4721875 ),
              ( 0.0946875 ,  0.5590625 ),
              ( 0.1415625 ,  0.6459375 ),
              ( 0.231875  ,  0.69234375),
              ( 0.30078125,  0.72796875),
              ( 0.40328125,  0.72796875),
              ( 0.5365625 ,  0.72796875),
              ( 0.6115625 ,  0.67203125),
              ( 0.6865625 ,  0.61625   ),
              ( 0.70796875,  0.51765625),
              ( 0.56453125,  0.49078125),
              ( 0.549375  ,  0.5434375 ),
              ( 0.5075    ,  0.57390625),
              ( 0.46578125,  0.60453125),
              ( 0.40328125,  0.60453125),
              ( 0.30859375,  0.60453125),
              ( 0.25265625,  0.544375  ),
              ( 0.19671875,  0.484375  ),
              ( 0.19671875,  0.36625   ),
              ( 0.19671875,  0.23875   ),
              ( 0.25328125,  0.175     ),
              ( 0.31      ,  0.11140625),
              ( 0.401875  ,  0.11140625),
              ( 0.44734375,  0.11140625),
              ( 0.49296875,  0.12921875),
              ( 0.53859375,  0.14703125),
              ( 0.57125   ,  0.17234375),
              ( 0.57125   ,  0.263125  )];

       t =   [(0.23390625, 0.        ),
              (0.23390625, 0.5946875 ),
              (0.0215625 , 0.5946875 ),
              (0.0215625 , 0.71578125),
              (0.5903125 , 0.71578125),
              (0.5903125 , 0.5946875 ),
              (0.3784375 , 0.5946875 ),
              (0.3784375 , 0.        )]

       _agct = Dict(DNA_A => a, DNA_C => [c], DNA_G => [g], DNA_T => [t])
       agct = Dict(b => normletter(l) for (b, l) in _agct)

       agct
end