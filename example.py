# coding: utf-8
import B4x4
sim = B4x4.Factory()
sim.points = 401
sim.nDelta = 0
result = sim.calculateL()
result.plot()
result.save('left_reflection.txt')

