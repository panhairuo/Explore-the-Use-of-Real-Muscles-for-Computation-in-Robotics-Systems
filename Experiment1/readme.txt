Description of Files

Here follows a short description of the folders and files:

• folder data_NARMA:
	– NARMA-L2.mat: example learning data
	– make_NARMA_data.m: Matlab code to produce learning databased on a NARMA (nonlinear autoregressive–moving-average) system


• folder helping_functions:
	– e_distance.m: calculates the euclidean distance between two points (to get current spring lengths)
	– learn_linear_model.m: learns a linear model with linear regression to emulate a given input/output data set (to evaluate the computational power of the 		 	   morphology）
	– make_simulation_movie.m: function to produce an .avi movie of the simulated network
	– nth_point.m: function used for plot_one_period.m (reduces the size of the figure file)
	– plot_graph.m: plots the networks with color coding (as in Fig. 7 of [1])
	– rand_in_range.m: draws random numbers from a given range
	– rand_in_range_exp.m: draws random numbers from a given range with an exponential distribution.
	– rhs_mass_sys.m: function used for the Runge-Kutta algorithm
	– rk4ode2.m: Runge-Kutta algorithm used for the physical simulation
	– step_response.m: produces a step response (or impulse response) of a given network
• main folder:
	– init_ms_sys_data.m: initializes a data structure with default values that defines the properties and parameter ranges for a network. The returned data structure can be 	   manipulated and is fed to init_ms_sys_net.m.
	– ode_simple_ms_sys.m: Implements a simple Euler integration scheme for a 2D mass-point system. Updates position and velocity states based on applied forces.

	– init_ms_sys_net.m: Initializes a mass-spring network structure with random parameters. Defines point masses, spring connections, and input/output weights. Includes 	   optional Bouc-Wen hysteresis parameters for advanced spring modeling.

	– plot_graph.m: Visualizes the mass-spring network structure. Plots mass points, spring connections, and highlights fixed, input, and output nodes.

	Three variants of the mass-spring network simulation function:
	– simulate_ms_sys_1.m: Basic nonlinear spring model (cubic + linear)
   	– simulate_ms_sys_2.m: Modified nonlinear model with max(0,x) activation
 	– simulate_ms_sys_3.m: Advanced model with Bouc-Wen hysteresis

	– Experiment1.m: Main experiment script that compares the three simulation variants. Loads NARMA dataset, trains networks, and evaluates performance. Measures 	   MSE and computation time for each model.

