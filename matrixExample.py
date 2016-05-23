# coding: utf-8
import B4x4
sim = B4x4.Factory()
result = sim.matrix('angle',0,90)
result.image()
result.save('matrix.txt')