# coding: utf-8
import B4x4
sim = B4x4.Simulator()
result_left = sim.calculateL()
result_left.plot()
result_right = sim.calculateR()
result_right.plot()
result_left.save('left_reflection.txt')
result_right.save('right_reflection.txt')
