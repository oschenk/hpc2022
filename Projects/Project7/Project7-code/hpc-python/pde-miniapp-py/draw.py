from pde_miniapp_py import viz

# initial
viz.draw("simulation", 0)
viz.plt.show()
viz.plt.savefig("simulation_initial.png")

# final
viz.draw("simulation", 100)
viz.plt.show()
viz.plt.savefig("simulation_T.png")

