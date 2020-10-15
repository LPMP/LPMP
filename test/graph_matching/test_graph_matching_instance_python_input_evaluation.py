import numpy as np
import graph_matching_py as gm

instance = gm.graph_matching_input();
# TODO: add instance to repo
instance.read("/local/home/pswoboda/test_LPMP/LPMP/test/graph_matching/instance_solverweirdness.txt")

solver = gm.graph_matching_message_passing_solver(['','--graphMatchingRounding','fw','--maxIter','3','--standardReparametrization','uniform:0.5'])
solver.construct(instance)
solver.solve()
result = solver.result()
print("instance evaluation after optimization: " + str(instance.evaluate(result)))
repam_instance = solver.export()
print("repam instance evaluation after optimization: " + str(repam_instance.evaluate(result)))
