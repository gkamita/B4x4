# coding: utf-8
import B4x4
sim = B4x4.Factory()
sim.points = 50
sim.lbda_max = 650
result = sim.matrix('angle',0,90)
result.image()
result.save('matrix.txt')
