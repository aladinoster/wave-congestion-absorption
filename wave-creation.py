"""
    Traffic Wave Absorption for CAVs

    @author Andres Ladino 

    Check README.md for more details
"""

import symupy
from symupy.api import Simulator
from symupy.components import VehicleControl
import numpy as np
import os


if __name__ == "__main__":

    # Modify the simulation insert points
    filename_path = os.path.join(os.getcwd(), "network/ring_48_cav.xml")
    sim_path = os.path.join(
        "/Users/ladino/Documents/03-Code/02-Python/libraries/symupy/lib/osx-64/libSymuVia.dylib"
    )
    sim_instance = Simulator.from_path(filename_path, sim_path)

    # Control
    c = VehicleControl()
    a = np.zeros(10)
    a[5] = -1
    a

    ac = iter(a)

    with sim_instance as s:
        while sim.do_next:
            sim.run_step()
            acc = next(ac)
            #         print(s.state.current_time)
            print(s.state.current_time == "5.00")
            print(s.state.vehicles)
            if sim.state.is_vehicle_in_network("0", "1") and sim.state.current_time == "7.00":
                c = VehicleControl("manual", 1, s.state.vehicles[0])
                c.set_manual_control(acc)
                sim.stop_step()


#     # 1. Declare V2V Network based on the vehicle type?

#     # 2. What is fundamental for the simulation is the XML file.
#     #    The idea is that the scenario should run independently of having modifycations via this API. so eventually you can just launch your simulation via this API.

#     # 3. Modify simulation such that she may accept a V2V object
#     #    If the simulation detects a network object the simulation
#     #    is enriched with the V2V characteristics.

#     # 4. Control means that
#     #    A control means a manipulation for 1 variable of a vehicle,       acceleration or speed

#     # Run simulation
#     simulation.run_Simulation()


# # """
# #     Fuel consumption analysis for Truck platooning

# #     This script runs an scenario of 4 trucks in platoon mood where
# #     particular splits are determined.

# #     In order to use: python platoon-closed.py

# #     Check output files in: ../Output/
# # """
# # from parameters import VehParameter, SimParameter, CtrParameter
# # from models import VehNetwork, Vehicle, dynamic_3rd, dynamic_2nd
# # from control import OperationalCtr, TacticalCtrl


# # # Length of the platoon
# # N_VEH = 3

# # # Create a simulation timings
# # T_STP = 0.01
# # T_HOR = 0.5
# # T_SIM = 60

# # sim_par = SimParameter(T_STP, T_HOR, T_SIM)

# # # Create a vehicle model / Provided originally by the simulator
# # U_FFS = 25.0
# # K_X = 0.16
# # W_CGT = 6.25
# # L_VEH = 4.0

# # S0 = 10.0

# # veh_par = VehParameter.VehParameterSym(U_FFS, K_X, W_CGT, L_VEH)

# # # Create list of vehicles
# # list_id = range(N_VEH)
# # veh_list = [Vehicle(sim_par, veh_par, dynamic_3rd, id=i) for i in list_id]

# # # Artificial leaders:
# # # Created a-priori
# # list_id = list(list_id)[1:]
# # lead_id = [i-1 for i in list_id[0:-1]]
# # net_veh = dict(zip(list_id, lead_id))

# # # Create the network of vehicles
# # platoon = VehNetwork(sim_par, veh_list)
# # print(f'Vehicle network: {net_veh}')
# # platoon.register_vehicle_link(net_veh)


# # # Initialize the network of vehicles
# # x0 = [0.0] * N_VEH
# # s0 = [S0] * N_VEH
# # v0 = [U_FFS] * N_VEH
# # e0 = [0.0] * N_VEH
# # a0 = [0.0] * N_VEH  # only 3rd order models

# # # Initialize each vehicle
# # state0 = [s0, v0, e0, a0]
# # state0veh = list(zip(*state0))  # state each vehicle
# # state0net = dict(zip(platoon.veh_currentids, state0veh))
# # platoon.initialize_vehicles(state0net)

# # # Create the controller[]
# # ctr_par = CtrParameter()

# # split_events = {1: {'ta': 20, 'tm': 40, 'tau0': 2, 'tauf': 5}, }
# # tc_ctrl = TacticalCtrl(platoon.sim_par, ctr_par, split_events)
# # tc_ctrl.compute_reference()

# # op_ctrl = OperationalCtr(platoon.sim_par, ctr_par)
# # op_ctrl.register_veh_network(platoon)
# # # op_ctrl.build_global_dynamics()
# # # op_ctrl.build_lagrange_dynamics()
# # op_ctrl.evolve_forward_dynamics()
# # # veh_network.launch_simulation()

# # Scenario 1: Perturb the measurements with noise

# # Scenario 2: Perturb the model parameters with noise / Parameters should keep consistancy

# # Scenario 3: Save on each vehicle the full buffer but remember to retrieve the info in the past. Compute ctrl with this info and then feed the system.

# # Obtain Performance measurements

# # In all cases: Dynamic Answer (State + Control)

# # Save results
