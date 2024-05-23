
import SuctionSkirtSolver.SuctionSkirtSolver as suskso
from plotting_files import main_plotter

if __name__=='__main__':
	setup_vars = suskso.load_yaml("setup_variables.yaml")
	skirt_stiffness_key = "stiffness_skirt"
	skirt_stiffness = [1,0.7,0.4,0.1]
	resolution = 10

	suskso.Full_Orientation_Test(resolution,setup_vars,override_key=skirt_stiffness_key,override_vars=skirt_stiffness,write=True)
	main_plotter()