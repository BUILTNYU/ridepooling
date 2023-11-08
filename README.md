## Introduction
NOMAD-RPS (Non-myopic Optimization for Matching and Allocation in Dynamic Ride-Pooling Simulation) is a dynamic ride-pooling simulator that uses a non-myopic cost function approximation policy to make sequential decisions and find the best matches between riders and vehicles in the system. Ride-pooling services promise to be more efficient than private motorized mobility or ride-hailing. It can reduce costs for both passengers and operators and reduce traffic congestion and environmental impacts. The static pickup and delivery problem with transfers is an NP-hard problem, while efficient heuristics are necessary for online applications due to the "curse of dimensionality". Myopic decision-making is another well-known issue in dynamic routing, which is often ignored. It is shown that there are opportunity costs to decisions made in routing at a point in time that can impact outcomes in the future horizon. Therefore, it is crucial to adopt a look-ahead approach in fleet operations. This simulator uses an online non-myopic policy and algorithm for operating a ride-pooling service which controls for opportunity costs due to commitment to serving existing passengers. 

## Problem description and methodology
The road network is a directed graph of N nodes and E edges, G = (N,E). The weight of each edge is the travel cost of that link in the network, which is represented by travel time. A set of nodes is defined as virtual stops where passengers can be picked up or dropped off. Requests are submitted during the service operating hours within the service region. Each trip request has the following information: submission time, origin location, and destination location. It is assumed that each request is associated with one passenger, and when a trip request is submitted, the passenger is ready for immediate pickup. Once the trip request is submitted, the origin and destination locations are assigned to nearby virtual stops. The initial location of vehicles can be randomly selected from the network nodes or initiated from predetermined hub locations. The system is updated after each time window. Requests submitted during a time window are sent to the centralized system for vehicle assignment. Sequential dispatching decisions are made with the goal of maximizing or minimizing an objective function that defines a performance criteria, such as passengers’ journey time and vehicles’ travel time. The following constraints are considered for finding a feasible match between a vehicle and a rider:
- Wait time for pickup from the origin for each passenger should be less than a maximum wait time.
- Drop-off time at the destination for each passenger should be less than their latest drop-off time.
- Vehicle capacity should not be violated at any point in time.
- Pickup time for each passenger should be earlier than their drop-off time.

The one-dimensional non-myopic policy for making decisions at each time step can be written as follows:
![image](https://github.com/BUILTNYU/ridepooling/assets/66441622/0868ebf7-9882-455e-929a-865d70a5d0a3)

where c(v, ξ) is the cost value of the routing/dispatching decision (vehicle v, route ξ), T(v, ξ) represents the operator cost, which is vehicle v’s travel time given route ξ. Parameter θ is the degree of operator cost vs. user cost. User cost is represented by passengers’ perceived travel time, which is a weighted sum of their in-vehicle time and wait time. J<sub>in−vehicle,n</sub>(v, ξ) is the in-vehicle travel time for passenger n assigned to vehicle v given route ξ, which is their drop-off time subtracted by their pickup time. J<sub>wait,n</sub>(v, ξ) is the wait time for passenger n assigned to vehicle v given route ξ, which is determined by subtracting the submission time from the pickup time. According to [Yap and Cats (2023)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4355442), waiting time for ride-hailing services is perceived about 1.4 times more negatively than the in-vehicle time. Therefore, a coefficient of 1.4 is chosen for wait time. β is the degree of lookahead approximation, which should be tuned for the problem.

The simulator also supports rebalancing of idle vehicles. This feature can be enabled in the code.

## Instructions
Simulation parameters and control variables can be determined in the code. Different input and output files are explained below:

### Input files
Example input files are available for the Sioux Falls network. The input files include:

- Network links, which contains the travel time on each link.
- Network nodes, which contains all the nodes in the network and their projected coordinates.
- Service stops, which contains the information of the nodes where passengers can be picked up or dropped off.
- Travel time matrix, which contains precalculated shortest travel time from each service stop to another.
- Hub locations, which contains the location of hubs from which vehicles start their operations (optional; the simulator can select the initial location of vehicles randomly from the network nodes).
- Requests, which contains the information of each submitted request including origin, destination, and submission time (optional; the simulator can generate random requests given the total number of submitted requests during operating hours).

### Output files
The output files include:

- Performance evaluation, which contains the values of performance metrics, including average passengers’ in-vehicle time, wait time, perceived journey time (includes in-vehicle time and wait time multiplied by a factor of 1.4), total vehicles’ travel time excluding the dwell times, total vehicles’ empty travel time which refers to the duration that vehicles are moving while not carrying any passengers, average vehicle occupancy, number of rejected requests, simulation runtime.
- Vehicle occupancy information at each time step.
- Requests' information, which contains each request's submission time, pickup time, drop-off time
- Vehicles' real-time location, which contains the location coordinates of every vehicle in the system at each time step.

## License
The NYU non-commercial research license is applied to NOMAD-RPS (attached in the repository). Please contact Joseph Chow ([joseph.chow@nyu.edu](joseph.chow@nyu.edu)) for commercial use.
For questions about the code, please contact Farnoosh Namdarpour ([farnoosh@nyu.edu](farnoosh@nyu.edu)).
