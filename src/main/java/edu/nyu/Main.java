package edu.nyu;

import com.opencsv.CSVWriter;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

public class Main {

    static final int big_N = (int) Math.pow(10,7);
    public static void main(String[] args) throws CloneNotSupportedException, IOException {

        //FILENAMES
        String filename = "SiouxFalls"; //this name is used for generating outputs
        String ttMatrix_filename = "SiouxFalls_TTmatrix.csv";
        String nodes_filename = "SiouxFalls_nodes.csv";
        String links_filename = "SiouxFalls_links.csv";
        String requests_filename = "SiouxFalls_requests.csv";
        String stops_filename = "SiouxFalls_stops.csv";
        String hubs_filename = "SiouxFalls_hubs.csv";

        //CONTROLS
        boolean use_hubs = true; //if true, the vehicles are initiated from hub locations otherwise they are assigned random initial locations
        boolean create_requests = false; //if true, it ignores requests_filename and generates random requests
        boolean write_occupancy_info = false;
        boolean write_requests_info = true;
        boolean enable_rebalancing = true;
        boolean write_realtime_vehicle_locations = true;

        //PARAMETERS
        int simulation_period = 24 * 3600; //seconds
        int simulation_start_time = 0 * 3600;
        int t_step = 30; //time step in terms of seconds
        int n_veh = 20; //number of vehicles
        int v_cap = 6; //vehicle capacity
        int n_nodes = 24; //number of nodes in the network
        int n_requests = 3000;
        int max_t_dist = 10 * 60; //maximum travel distance in terms of time unit (seconds) for finding candidate vehicles
        int max_wait_t = 10 * 60; //maximum waiting time for passengers to be picked up
        float max_tt_p = 0.4f; //a parameter for finding maximum travel time for each request, provided by MOIA
        int max_tt_min = 10 * 60; //a parameter for finding maximum travel time for each request, provided by MOIA
        int max_tt_added = 15 * 60; //a parameter for finding maximum travel time for each request, provided by MOIA - maximum travel time added to the direct travel time
        int base_dwell_time = 50; //base dwell time at each stop (seconds)
        int p_dwell_time = 10; //dwell time added to base dwell time for each passenger picked up or dropped off (seconds)
        float wait_vot = 1.4f; //value of waiting time relative to in-vehicle travel time
        float theta = 0.5f; //a parameter of the cost function
        float beta = 0.005f; //a parameter of the cost function
        double big_M = Double.POSITIVE_INFINITY; //an arbitrary large value used in different functions
//        long seed = 4; //random seed
        String requests_source_crs = "EPSG:4326"; //EPSG:4326 is equivalent to WGS84
        String requests_destination_crs = "EPSG:3455"; //EPSG:25832 is projected coordinate system in meters for Europe, EPSG:3455 is in ft for SiouxFalls
        String hubs_source_crs = "EPSG:4326"; //EPSG:4326 is equivalent to WGS84
        String hubs_destination_crs = "EPSG:3455"; //EPSG:25832 is projected coordinate system in meters for Europe, EPSG:3455 is in ft for SiouxFalls
        String nodes_source_crs = "EPSG:4326"; //EPSG:4326 is equivalent to WGS84
        String nodes_destination_crs = "EPSG:3455"; //EPSG:25832 is projected coordinate system in meters for Europe, EPSG:3455 is in ft for SiouxFalls


        HashMap<Integer, Vehicle> vehicles = new HashMap<>();
        HashMap<Integer, Request> requests = new HashMap<>();
        Data data = new Data(new HashMap<>()); //used for storing simulation data
        Network network = new Network(n_nodes);
        HashMap<Integer, Node> nodes = new HashMap<>();
        HashMap<Integer, Node> stops = new HashMap<>();
        ArrayList<Float> avg_occ = new ArrayList<>(); //vehicles' average occupancy
        ArrayList<Integer> stop_ids = new ArrayList<>();
        ArrayList<Integer> hub_locations = new ArrayList<>();

        //write performance evaluation results
        File file = new File(filename + "_results.csv");
        FileWriter outputfile = new FileWriter(file);
        CSVWriter writer = new CSVWriter(outputfile);
        String[] header = {"seed", "simulation_period", "n_requests", "n_veh", "n_stops", "theta",
                "beta", "avg pax in-vehicle time (min)", "avg pax waiting time (min)",
                "avg pax travel time (min)", "vehicles idle time (hrs)", "avg vehicle occupancy", "# of pax served",
                "# of rejected or pending", "total traveled distance (hrs)", "avg traveled distance (min)",
                "total empty traveled distance (hrs)", "run time"};
        writer.writeNext(header);

                //data for occupancy plot
                File file_occ = new File(filename + "_occupancyInfo.csv");
                FileWriter outputfile_occ = new FileWriter(file_occ);
                CSVWriter writer_occ = new CSVWriter(outputfile_occ);
                String[] header_occ = {"time", "idle", "occ_0", "occ_1", "occ_2", "occ_3", "occ_4", "occ_5", "occ_6"};
                writer_occ.writeNext(header_occ);

                //write realtime vehicle locations
                File file_vl = new File(filename + "_v_realtime_locations.csv");
                FileWriter outputfile_vl = new FileWriter(file_vl);
                CSVWriter writer_vl = new CSVWriter(outputfile_vl);
                String[] header_vl = new String[n_veh*2+1];
                header_vl[0]="time(sec)";
                for (int g=1; g<header_vl.length; g=g+2){
                    header_vl[g]="V" + (g-1)/2 + "_X";
                    header_vl[g+1] = "V" + (g-1)/2 + "_Y";
                }
                writer_vl.writeNext(header_vl);

                for (int seed = 0; seed < 1; seed++) {

                    long start_time = System.nanoTime();

                    try {
                        System.out.println("reading network");
                        Reader.read_network(network, ttMatrix_filename);
                        System.out.println("reading nodes");
                        Reader.read_nodes(nodes, nodes_filename, nodes_source_crs, nodes_destination_crs);
                        System.out.println("reading stops");
                        Reader.read_stops(stop_ids, nodes, stops, stops_filename, nodes_source_crs, nodes_destination_crs);
                        System.out.println("reading requests");
                        if (create_requests) {
                            create_requests(n_requests, requests, stop_ids, simulation_period, network, max_wait_t, max_tt_p, max_tt_min, max_tt_added, seed);
                        } else {
                            Reader.read_requests(requests, stops, network, requests_filename, max_wait_t,
                                    max_tt_p, max_tt_min, max_tt_added, requests_source_crs, requests_destination_crs);
                        }
                        n_requests = requests.size();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }

                    //create a graph for shortest path calculation
                    DirectedGraph graph = new DirectedGraph();
                    create_graph(links_filename, n_nodes, graph);

                    System.out.println("number of nodes in the network: " + nodes.size());
                    System.out.println("number of stops in the network: " + stop_ids.size());
                    System.out.println("number of requests: " + requests.size());

                    // INITIALIZE VEHICLE LOCATIONS
                    System.out.println("initializing vehicles location");
                    if (use_hubs) {
                        Reader.read_hub_locations(hub_locations, stops, hubs_filename,
                                hubs_source_crs, hubs_destination_crs);
                    }
                    System.out.println("hub locations size: " + hub_locations.size());

                    initialize_vehicles(use_hubs, n_veh, stop_ids, vehicles, seed, hub_locations);
                    System.out.println("vehicles size:" + vehicles.size());

                    // INITIALIZE DATA FOR VEHICLES
                    System.out.println("initializing data for vehicles");
                    initialize_sim_data(vehicles, data);

                    //SEPARATE REQUESTS BY TIME STEPS
                    System.out.println("separating requests by time steps");
                    ArrayList<Integer>[] ordered_requests;
                    ordered_requests = separate_requests(requests, t_step, simulation_period, simulation_start_time);

                    //PROCESS REQUESTS IN EACH TIME STEP
                    int n = simulation_period / t_step; //number of time steps in the simulation period
                    for (int x = 0; x < n; x++) {
                        int t = x * t_step + simulation_start_time;
                        System.out.println("t: " + t);

                        if (write_occupancy_info) {
                            int occ[] = new int[]{0, 0, 0, 0, 0, 0, 0, 0};
                            for (Map.Entry<Integer, Vehicle> vehicle : vehicles.entrySet()) {
                                if (vehicle.getValue().route_nodes.size() == 0 && vehicle.getValue().tt_to_next_node == 0) {
                                    occ[0] += 1;
                                } else if (vehicle.getValue().current_load == 0) {
                                    occ[1] += 1;
                                } else if (vehicle.getValue().current_load == 1) {
                                    occ[2] += 1;
                                } else if (vehicle.getValue().current_load == 2) {
                                    occ[3] += 1;
                                } else if (vehicle.getValue().current_load == 3) {
                                    occ[4] += 1;
                                } else if (vehicle.getValue().current_load == 4) {
                                    occ[5] += 1;
                                } else if (vehicle.getValue().current_load == 5) {
                                    occ[6] += 1;
                                } else if (vehicle.getValue().current_load == 6) {
                                    occ[7] += 1;
                                }
                            }
                            String[] line_occ = {String.valueOf(x), String.valueOf(occ[0]), String.valueOf(occ[1]), String.valueOf(occ[2]),
                                    String.valueOf(occ[3]), String.valueOf(occ[4]), String.valueOf(occ[5]), String.valueOf(occ[6]), String.valueOf(occ[7])};

                            writer_occ.writeNext(line_occ);
                        }

                        time_update(theta, beta, t, t_step, vehicles, data, requests, network, base_dwell_time, p_dwell_time, max_wait_t, avg_occ, wait_vot);

                        for (int r_id : ordered_requests[x]) {

//                print_request_info(requests, r_id);

                            //find candidate vehicles for finding either best vehicle without transfers or best first vehicle for transfer
                            ArrayList<Integer> candidate_vehicles = find_candidate_vehicles(requests.get(r_id), max_t_dist, vehicles, network);

                            //find the best vehicle without transfer
                            Container best_vehicle; //container format: veh_id, cost, route

                            best_vehicle = find_best_vehicle(network, vehicles, requests, candidate_vehicles, t, r_id, v_cap,
                                    max_wait_t, theta, beta, p_dwell_time, base_dwell_time, big_M, t_step, wait_vot);

                            Integer best_v_id = best_vehicle.v1_id;
                            List<Vehicle.Route_stop> best_v_route = best_vehicle.v1_route;
                            Double best_v_cost = best_vehicle.v1_cost;
                            double best_cost_diff = big_M;
                            if (best_v_id != null) {
                                best_cost_diff = best_v_cost - vehicles.get(best_v_id).cost;
                            }

//                print_transfer_info(r_id, 549, requests, vehicles, best_v_id, best_v_route, best_cost_diff, best_transfer, transfer_cost);

                            //handling requests that had no feasible solution
//                            if (best_vehicle.v1_id == null && best_transfer.v2_id == null) {
                            if (best_vehicle.v1_id == null) {
                                //TODO should it be t or t+t_step?
                                if (t - requests.get(r_id).submission_t > max_wait_t) {
                                    requests.get(r_id).status = "rejected";
                                    if (enable_rebalancing) {
                                        //find the best vehicle to go to the pickup location of rejected request
                                        Container best_vehicle_rebalance; //container format: veh_id, cost, route

                                        best_vehicle_rebalance = find_best_vehicle_rebalance(network, vehicles, requests, t, r_id, v_cap,
                                                max_wait_t, theta, beta, p_dwell_time, base_dwell_time, big_M, t_step, wait_vot);

                                        if (best_vehicle_rebalance.v1_id != null) {
                                            Integer best_v_id_rebalance = best_vehicle_rebalance.v1_id;
                                            List<Vehicle.Route_stop> best_v_route_rebalance = best_vehicle_rebalance.v1_route;
                                            Double best_v_cost_rebalance = best_vehicle_rebalance.v1_cost;
                                            assignment_update(vehicles, requests, network, best_v_id_rebalance, best_v_route_rebalance, r_id, best_v_cost_rebalance, graph, true);
                                        }
                                    }
                                } else {
                                    requests.get(r_id).status = "pending";
                                    if (x + 1 < ordered_requests.length) {
                                        ordered_requests[x + 1].add(r_id);
                                    }
                                }
                            }  else {
                                assignment_update(vehicles, requests, network, best_v_id, best_v_route, r_id, best_v_cost, graph, false);
                            }

//                print_route(vehicles, 9, new ArrayList<>());
//                System.out.println("cost: " + vehicles.get(9).cost);

                        }


//            System.out.println("after assignment");
//            print_route(vehicles, 9, new ArrayList<>());
//            System.out.println("cost: " + vehicles.get(9).cost);

                        if (write_realtime_vehicle_locations) {
                            String[] line_vl = new String[n_veh*2+1];
                            line_vl[0] = Integer.toString(t);
                            for (int g=1; g<line_vl.length; g=g+2){
                                line_vl[g] = Double.toString(nodes.get(vehicles.get((g-1)/2).current_loc).x);
                                line_vl[g+1] = Double.toString(nodes.get(vehicles.get((g-1)/2).current_loc).y);
                            }
                            writer_vl.writeNext(line_vl);
                        }
                    }

                    //one last time update (t is equal to simulation_period here)
                    time_update(theta, beta, simulation_period, t_step, vehicles, data, requests, network, base_dwell_time, p_dwell_time, max_wait_t, avg_occ, wait_vot);

                    //reporting simulation outputs
                    report_outputs(requests, data, simulation_period, n_veh, avg_occ, n, theta, beta, n_requests, stops.size(), start_time, writer, seed, wait_vot);

//
////                    print_transfers(requests, data);

//        print_request_info(requests, 170);
//        print_request_info(requests, 10170);


                    writer.flush();

                    //requests info output file
                    if (write_requests_info) {
                        report_requests_info(filename + "_requests_info.csv", requests, nodes);
                    }
                }
                writer_occ.close();
                if (write_occupancy_info == false) {
                    file_occ.delete();
                }
                writer_vl.close();
                if (write_realtime_vehicle_locations == false) {
                    file_vl.delete();
                }
        writer.close();
    }

    public static void create_graph(String links_filename, int n_nodes, DirectedGraph graph) {

        for (int i = 0; i < n_nodes; i++) {
            graph.addNode(i);
        }

        //Read Link info
        BufferedReader br;
        try {
            br = new BufferedReader(new FileReader(links_filename));
            String line;
            br.readLine();
            while ((line = br.readLine()) != null) {
                String[] nn_data = line.split(",");

                int src_node = (int) Double.parseDouble(nn_data[10]);
                int dest_node = (int) Double.parseDouble(nn_data[11]);
                double duration = Double.parseDouble(nn_data[12]);
                graph.addEdge(src_node, dest_node, duration);
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static List<Integer> find_shortest_path(int origin, int dest, DirectedGraph graph) {

        Map<Integer, Integer> path = new HashMap<>();
        Map<Integer, Double> result = new HashMap<>();
        HashMap<Map<Integer, Double>, Map<Integer, Integer>> output;

        output = Dijkstra.shortestPaths(graph, origin, dest);

        for (Map.Entry<Map<Integer, Double>, Map<Integer, Integer>> pair : output.entrySet()) {
            path = pair.getValue();
            result = pair.getKey();
        }

        int a = dest;
        List<Integer> shortestPath = new ArrayList<>();

        while (a != origin) {
            shortestPath.add(a);
            a = path.get(a);
        }
        shortestPath.add(origin);

        Collections.reverse(shortestPath);
        return shortestPath;
    }

    public static void initialize_vehicles(boolean use_hubs, int n_veh, List<Integer> stop_ids, HashMap<Integer, Vehicle> vehicles, long seed, ArrayList<Integer> hub_locations) {

        if (use_hubs) {
            for (int i = 0; i < hub_locations.size(); i++) {
                for (int j = i * n_veh / hub_locations.size(); j < (i + 1) * n_veh / hub_locations.size(); j++) {
                    int loc = hub_locations.get(i);
                    Vehicle vehicle = new Vehicle(loc, 0, 0, 0, new ArrayList<>(), new ArrayList<>(), new ArrayList<>());
                    vehicles.put(j, vehicle);
//                    System.out.println("j: " + j);
//                    System.out.println("location: " + loc);
                }
            }
        } else {
            Random rand = new Random(seed);
            for (int i = 0; i < n_veh; i++) {
                int loc = stop_ids.get(rand.nextInt(stop_ids.size()));
                Vehicle vehicle = new Vehicle(loc, 0, 0, 0, new ArrayList<>(), new ArrayList<>(), new ArrayList<>());
                vehicles.put(i, vehicle);
            }
        }
    }

    public static void initialize_sim_data(HashMap<Integer, Vehicle> vehicles, Data data) {

        for (int key : vehicles.keySet()) {
            data.vehicles.put(key, new Data.VehicleSim(0, 0, 0, new ArrayList<>(), 0, vehicles.get(key).current_loc, 0));
        }
    }

    public static void create_requests(int n_requests, HashMap<Integer, Request> requests, List<Integer> stop_ids,
                                       int simulation_period, Network network, int max_wait_t, float max_tt_p, int max_tt_min,
                                       int max_tt_added, long seed) {

        Random rand = new Random(seed);
        for (int i = 0; i < n_requests; i++) {
            int origin = stop_ids.get(rand.nextInt(stop_ids.size()));
            int dest = stop_ids.get(rand.nextInt(stop_ids.size()));
            while (origin == dest)
                dest = stop_ids.get(rand.nextInt(stop_ids.size()));
            int sub_time = rand.nextInt(simulation_period-1*3600);
            requests.put(i, load_request(origin, dest, sub_time, network.dist, max_wait_t, max_tt_p, max_tt_min, max_tt_added));
        }
    }

    public static Request load_request(int origin, int dest, int submission_t, int[][] dist, int max_wait_t,
                                       float max_tt_p, int max_tt_min, int max_tt_added) {

        int direct_tt = dist[origin][dest];
        int latest_arr_t;
        if (max_tt_p * direct_tt + max_tt_min > max_tt_added) {
            latest_arr_t = (int) (submission_t + max_wait_t + direct_tt + max_tt_added);
        } else {
            latest_arr_t = (int) (submission_t + max_wait_t + (1 + max_tt_p) * direct_tt + max_tt_min);
        }
        Request request = new Request(origin, dest, submission_t, latest_arr_t, null, null,
                "submitted", new ArrayList<>(), false);
        return request;
    }

    public static ArrayList<Integer>[] separate_requests(HashMap<Integer, Request> requests, int t_step, int simulation_period, int simulation_start_time) {

        int n = simulation_period / t_step;
        ArrayList<Integer>[] ordered_requests = new ArrayList[n];

        // initializing
        for (int i = 0; i < n; i++) {
            ordered_requests[i] = new ArrayList<>();
        }

        for (Map.Entry<Integer, Request> request : requests.entrySet()) {

            for (int t = simulation_start_time; t < simulation_period + simulation_start_time; t += t_step) {
                if (request.getValue().submission_t >= t & request.getValue().submission_t < t + t_step) {
                    ordered_requests[(t - simulation_start_time) / t_step].add(request.getKey());
                    break;
                }
            }
        }

        return ordered_requests;
    }

    public static void time_update(float theta, float beta, int t, int t_step, HashMap<Integer, Vehicle> vehicles, Data data,
                                  HashMap<Integer, Request> requests, Network network, int base_dwell_time, int p_dwell_time,
                                  int max_wait_t, ArrayList<Float> avg_occ, float wait_vot) throws CloneNotSupportedException {

        int sum_occ = 0;

        for (Map.Entry<Integer, Vehicle> pair : vehicles.entrySet()) {

            //when a new request is assigned to a vehicle with empty route_nodes, the first node in route_nodes is the
            // current location which should be removed if it's not a stop on the route
            //TODO it would be less costly to have this in assignment_update...
            if (pair.getValue().route_nodes.size() > 0 &&
                    pair.getValue().current_loc == pair.getValue().route_nodes.get(0) &&
                    !pair.getValue().route_nodes.get(0).equals(pair.getValue().route_stops.get(0).stop_id)) {
                pair.getValue().route_nodes.remove(0);
                //tt_to_next_node shows travel time to next node at the end of a time step. So if t=100, it shows tt_to_next_node at
                //the end of t+t_step which would be t=130
                pair.getValue().tt_to_next_node += network.dist[pair.getValue().current_loc][pair.getValue().route_nodes.get(0)];
            }

            //update idle time if vehicle has no other nodes to visit and doesn't need to wait at its current location
            if (pair.getValue().tt_to_next_node == 0 && pair.getValue().route_nodes.size() == 0) {
                data.vehicles.get(pair.getKey()).total_idle_t += t_step;

                //only update tt_to_next_node when vehicle doesn't reach next node in this time step
                //this condition also takes care of dwell/transfer time at final stops
            } else if (pair.getValue().tt_to_next_node > t_step) {
                pair.getValue().tt_to_next_node -= t_step;

                //in this case both current location and tt_to_next_node should be updated
            } else {
                int remaining_t_step = t_step - pair.getValue().tt_to_next_node;
                boolean cont_loop = true;
                //this loop is added to allow visiting more than one node in one time step (travel time from one node to another might be less than 30 seconds)
                while (cont_loop) {

                    //update current location
                    if (pair.getValue().route_nodes.size() > 0) {
                        pair.getValue().current_loc = pair.getValue().route_nodes.get(0);
                        //remove the node from route_nodes
                        pair.getValue().route_nodes.remove(0);
                    }

                    //update travel time to next node
                    int idle_time = 0;
                    if (pair.getValue().route_nodes.size() > 0) {
                        pair.getValue().tt_to_next_node = network.dist[pair.getValue().current_loc][pair.getValue().route_nodes.get(0)] -
                                (remaining_t_step);

                        //if there is no next node, update idle time and set travel time to next node to zero
                    } else {
                        idle_time = remaining_t_step;
                        pair.getValue().tt_to_next_node = 0;
                        cont_loop = false;
                    }

                    //update stop visits
                    if (pair.getValue().route_stops.size() > 0 && pair.getValue().current_loc == pair.getValue().route_stops.get(0).stop_id) {

                        //handle pickups
                        for (int i : pair.getValue().route_stops.get(0).pickup_ids) {
                            requests.get(i).status = "onboard";
                        }
                        //handle dropoffs
                        for (int i : pair.getValue().route_stops.get(0).dropoff_ids) {
                            requests.get(i).status = "alighted";
                            pair.getValue().active_requests.remove(Integer.valueOf(i));
                        }

                        //move the visited node to simulation data
                        data.vehicles.get(pair.getKey()).n_visited_stops += 1;
                        data.vehicles.get(pair.getKey()).n_served_pax += pair.getValue().route_stops.get(0).dropoff_ids.size();

                        //update total traveled distance
                        if (data.vehicles.get(pair.getKey()).n_visited_stops == 1) {
                            data.vehicles.get(pair.getKey()).total_distance += network.dist[data.vehicles.get(pair.getKey()).initial_loc][pair.getValue().route_stops.get(0).stop_id];

                            //update empty traveled distance
                            if (pair.getValue().current_load == 0) {
                                data.vehicles.get(pair.getKey()).total_empty_distance += network.dist[data.vehicles.get(pair.getKey()).initial_loc][pair.getValue().route_stops.get(0).stop_id];
                            }

                        } else {
                            data.vehicles.get(pair.getKey()).total_distance += network.dist[data.vehicles.get(pair.getKey()).route_stops.get(data.vehicles.get(pair.getKey()).route_stops.size() - 1).stop_id][pair.getValue().route_stops.get(0).stop_id];

                            //update empty traveled distance
                            if (pair.getValue().current_load == 0) {
                                data.vehicles.get(pair.getKey()).total_empty_distance += network.dist[data.vehicles.get(pair.getKey()).route_stops.get(data.vehicles.get(pair.getKey()).route_stops.size() - 1).stop_id][pair.getValue().route_stops.get(0).stop_id];
                            }
                        }
                        data.vehicles.get(pair.getKey()).route_stops.add(pair.getValue().route_stops.get(0).clone());

                        //change vehicle load
                        pair.getValue().current_load += pair.getValue().route_stops.get(0).pickup_ids.size() -
                                pair.getValue().route_stops.get(0).dropoff_ids.size();

                        //add dwell time to tt_to_next_node
                        //TODO if vehicle's current location is the same as stop id (maybe visited the stop and waiting due to dwell time) and another stop is added to the route with the same id, is it okay to reconsider base_dwell_time?
                        int dwell_time = base_dwell_time + (pair.getValue().route_stops.get(0).pickup_ids.size() + pair.getValue().route_stops.get(0).dropoff_ids.size()) * p_dwell_time;
                        pair.getValue().tt_to_next_node += dwell_time;

                        pair.getValue().tt_to_next_node -= idle_time;

                        //remove the stop from route_stops
                        pair.getValue().route_stops.remove(0);

                        //although, there might be no next nodes, vehicle might be waiting in the last stop due to remaining
                        // dwell/transfer time. In such cases, there is no idle time (base dwell time > time step)
                        idle_time = 0;
                    }
                    data.vehicles.get(pair.getKey()).total_idle_t += idle_time;
                    if (pair.getValue().tt_to_next_node < 0) {
                        remaining_t_step = -1 * pair.getValue().tt_to_next_node;
                    } else {
                        cont_loop = false;
                    }
                }
            }

            //updating cost
            HashMap<Integer, List<Integer>> active_requests_for_cost = new HashMap<>();
            for (int i : pair.getValue().active_requests) {
                active_requests_for_cost.put(i, Arrays.asList(requests.get(i).pickup_t, requests.get(i).dropoff_t, requests.get(i).submission_t));
            }
            pair.getValue().cost = find_cost(theta, beta, t + t_step, pair.getValue().route_stops, active_requests_for_cost,
                    base_dwell_time, p_dwell_time, wait_vot);

            //calculating average occupancy
            sum_occ += pair.getValue().current_load;
        }
        avg_occ.add((float) sum_occ / vehicles.size());
    }

    public static double find_cost(float theta, float beta, int t, List<Vehicle.Route_stop> route,
                                   HashMap<Integer, List<Integer>> v_requests, int base_dwell_time, int p_dwell_time, float wait_vot) {

        double cost = 0;
        if (route.size() > 0) {

            //find dwell time at the last stop for calculating vehicle's travel time
            int dwell_time = base_dwell_time + (route.get(route.size() - 1).pickup_ids.size() + route.get(route.size() - 1).dropoff_ids.size()) * p_dwell_time;

            int T = route.get(route.size() - 1).exp_stop_arr_t + dwell_time - t; //vehicle's travel time
            int S = 0; //sum of passengers' in-vehicle time
            int W = 0; //sum of passengers' waiting time
            for (Map.Entry<Integer, List<Integer>> pair : v_requests.entrySet()) {
                S += pair.getValue().get(1) - pair.getValue().get(0); //(0): pickup time, (1): dropoff time, (2): submission time
                W += pair.getValue().get(0) - pair.getValue().get(2);
            }
                cost = theta * T + (1 - theta) * S + wait_vot * (1 - theta) * W + (1 - theta) * beta * T * T;

        }
        return cost;
    }

    public static ArrayList<Integer> find_candidate_vehicles(Request request, int max_t_dist,
                                                             HashMap<Integer, Vehicle> vehicles, Network network) {

        ArrayList<Integer> candidate_vehicles = new ArrayList<>();
        //add all vehicles within a maximum travel time distance from the request's origin to the candidate_vehicles
        for (Map.Entry<Integer, Vehicle> pair : vehicles.entrySet()) {
            if (network.dist[pair.getValue().current_loc][request.origin] < max_t_dist) {
                candidate_vehicles.add(pair.getKey());
            }
        }
        return candidate_vehicles;
    }


    public static Container find_best_vehicle(Network network, HashMap<Integer, Vehicle> vehicles,
                                              HashMap<Integer, Request> requests, ArrayList<Integer> candidate_vehicles,
                                              int t, int r_id, int v_cap, int max_wait_t, float theta, float beta,
                                              int p_dwell_time, int base_dwell_time, double big_M, int t_step, float wait_vot) throws CloneNotSupportedException {

        double best_cost_diff = big_M;
        Container best_v = new Container(null, null, new ArrayList<>());
//        System.out.println("finding the best vehicle for request " + r_id);

        for (int v_id : candidate_vehicles) {
            HashMap<Double, List<Vehicle.Route_stop>> candidate_route; //the key is the cost of the route
            double cost_diff = big_M;
            double new_cost = big_M;
            List<Vehicle.Route_stop> route = new ArrayList<>();

            //insertion_heuristic finds the best insertion of pickup and dropoff nodes into the vehicle's route
            candidate_route = insertion_heuristic(t, v_id, vehicles.get(v_id), r_id, requests, network, v_cap, max_wait_t,
                    theta, beta, p_dwell_time, base_dwell_time, big_M, t_step, vehicles.size(), wait_vot);

            //candidate route should only have one item - hashmap is just used to return more than one object
            for (Map.Entry<Double, List<Vehicle.Route_stop>> item : candidate_route.entrySet()) {
                new_cost = item.getKey();
                cost_diff = new_cost - vehicles.get(v_id).cost;
                route = item.getValue();
            }
            if (cost_diff < best_cost_diff) {
                best_cost_diff = cost_diff;
                best_v = new Container(v_id, new_cost, route);
            }
        }

        //Merge repetitive consecutive stops
        merge_stops(best_v.v1_route);

        return best_v;
    }

    public static Container find_best_vehicle_rebalance(Network network, HashMap<Integer, Vehicle> vehicles,
                                                        HashMap<Integer, Request> requests, int t, int r_id, int v_cap,
                                                        int max_wait_t, float theta, float beta, int p_dwell_time,
                                                        int base_dwell_time, double big_M, int t_step, float wait_vot) throws CloneNotSupportedException {

        double best_cost_diff = big_M;
        Container best_v = new Container(null, null, new ArrayList<>());

        //find idle vehicles
        ArrayList<Integer> idle_vehicles = new ArrayList<>();
        for (Map.Entry<Integer, Vehicle> v : vehicles.entrySet()) {
            if (v.getValue().active_requests.size() == 0 && v.getValue().route_stops.size() == 0) {
                idle_vehicles.add(v.getKey());
            }
        }
        for (int v_id : idle_vehicles) {
            List<Vehicle.Route_stop> route = new ArrayList<>();

            //adding pickup node
            route.add(0, new Vehicle.Route_stop(requests.get(r_id).origin, null, new ArrayList<>(Arrays.asList()), new ArrayList<>(Arrays.asList())));

            update_exp_arr_t(0, base_dwell_time, p_dwell_time, t, t_step, network, vehicles.get(v_id), route, requests, max_wait_t);

            //temp_requests is used to have a copy of active requests and to make changes to them
            HashMap<Integer, List<Integer>> temp_requests = new HashMap<>();

            double cost = find_cost(theta, beta, t + t_step, route, temp_requests, base_dwell_time, p_dwell_time, wait_vot);

            //update the best route
            double cost_diff = cost - vehicles.get(v_id).cost;

            if (cost_diff < best_cost_diff) {
                best_cost_diff = cost_diff;
                best_v = new Container(v_id, cost, route);
            }
        }
        return best_v;
    }

    public static HashMap<Double, List<Vehicle.Route_stop>> insertion_heuristic(int t, int v_id, Vehicle v, int r_id,
                                                                                HashMap<Integer, Request> requests,
                                                                                Network network, int v_cap, int max_wait_t,
                                                                                float theta, float beta, int p_dwell_time,
                                                                                int base_dwell_time, double big_M,
                                                                                int t_step, int n_veh, float wait_vot) throws CloneNotSupportedException {
        double cost = big_M;
        int origin = requests.get(r_id).origin;
        int dest = requests.get(r_id).dest;
        HashMap<Double, List<Vehicle.Route_stop>> best_route = new HashMap<>(); //the key is the cost

        for (int i = 0; i <= v.route_stops.size(); i++) {
            for (int j = i + 1; j <= v.route_stops.size() + 1; j++) {

//                List<Vehicle.Route_stop> route = new ArrayList<Vehicle.Route_stop>(v.route_stops);
                List<Vehicle.Route_stop> route = new ArrayList<>();
                for (Vehicle.Route_stop item : v.route_stops) {
//                    route.add(Vehicle.Route_stop(item));
                    route.add(item.clone());
                }
                //add the pickup node
                route.add(i, new Vehicle.Route_stop(origin, null, new ArrayList<>(Arrays.asList(r_id)), new ArrayList<>(Arrays.asList())));
                //add the dropoff node
                route.add(j, new Vehicle.Route_stop(dest, null, new ArrayList<>(Arrays.asList()), new ArrayList<>(Arrays.asList(r_id))));


                merge_stops(route);

                //temp_requests is used to have a copy of active requests and to make changes to them
                HashMap<Integer, List<Integer>> temp_requests = new HashMap<>();

                //check feasibility
                int feasibility = check_feasibility(v, requests, r_id, route, temp_requests, v_cap, base_dwell_time,
                        p_dwell_time, t, t_step, network, max_wait_t);

                if (feasibility == 0) {
                    continue;
                } else {
                    double temp_cost = find_cost(theta, beta, t + t_step, route, temp_requests, base_dwell_time, p_dwell_time, wait_vot);

                    //update the best route
                    if (temp_cost < cost) {
                        cost = temp_cost;
                        best_route.clear();
                        List<Vehicle.Route_stop> route_copy = new ArrayList<>();
                        for (Vehicle.Route_stop item : route) {
                            route_copy.add(item.clone());
                        }
                        best_route.put(cost, route_copy);
                    }
                }
            }
        }
        return best_route;
    }

    public static void merge_stops(List<Vehicle.Route_stop> v_route) {

        for (int i = v_route.size() - 1; i > 0; i--) {
            int j = i - 1;
            if (v_route.get(i).stop_id.equals(v_route.get(j).stop_id)) {
                v_route.get(j).pickup_ids.addAll(v_route.get(i).pickup_ids);
                v_route.get(j).dropoff_ids.addAll(v_route.get(i).dropoff_ids);
                v_route.remove(i);
            }
        }
    }

    public static void update_exp_arr_t(int k, int base_dwell_time, int p_dwell_time, int t, int t_step, Network network,
                                        Vehicle v, List<Vehicle.Route_stop> route,
                                        HashMap<Integer, Request> requests, int max_wait_t) {

        int dwell_time = 0;
        if (k == 0) {
            //the second condition is added for the case when in a single time step, two requests are assigned to the
            // same vehicle which had no stops before. After the assignment of the first request to the vehicle, tt_to_next_node
            //is zero but route_nodes is not empty anymore.
            if (v.route_nodes.size() == 0 || v.tt_to_next_node == 0) {
                //v.tt_to_next_node is not necessarily zero in this case due to remaining dwell/transfer time at a final stop
                route.get(k).exp_stop_arr_t = t + t_step + v.tt_to_next_node + network.dist[v.current_loc][route.get(k).stop_id];
            } else {
                route.get(k).exp_stop_arr_t = t + t_step + v.tt_to_next_node + network.dist[v.route_nodes.get(0)][route.get(k).stop_id];
            }
        } else {

            if (k >= 2) {
                //in this case, the base_dwell_time has been considered in k-2th stop
                if (route.get(k - 2).stop_id.equals(route.get(k - 1).stop_id))
                    dwell_time = (route.get(k - 1).pickup_ids.size() + route.get(k - 1).dropoff_ids.size()) * p_dwell_time;
                else
                    dwell_time = base_dwell_time + (route.get(k - 1).pickup_ids.size() + route.get(k - 1).dropoff_ids.size()) * p_dwell_time;
            } else {
                dwell_time = base_dwell_time + (route.get(k - 1).pickup_ids.size() + route.get(k - 1).dropoff_ids.size()) * p_dwell_time;
            }
            //because of transfer waiting time considerations, dwell time is added later
            route.get(k).exp_stop_arr_t = route.get(k - 1).exp_stop_arr_t +
                    network.dist[route.get(k - 1).stop_id][route.get(k).stop_id];
        }

        //find maximum transfer waiting time in the previous stop
        //NOTE: if transfer stop is the final stop on the route, transfer waiting time is taken care of in the time_update (using tt_to_next_node)
        if (k >= 1) {
            int transfer_wait_t = 0;
            for (int p_id : route.get(k - 1).pickup_ids) {
                //when transfer is being evaluated in the find_best_transfer function, there is no need to add transfer_wait_t
                // It is taken care of in the insertion heuristic process
                // the first condition is for the case that transfer is being evaluated but not yet planned -> r_id+big_N is not in the requests object
                if (requests.containsKey(p_id) && requests.get(p_id).transfer == true && p_id > big_N) {
                    int temp_transfer_t = requests.get(p_id).submission_t - route.get(k - 1).exp_stop_arr_t + max_wait_t;
                    if (temp_transfer_t > transfer_wait_t)
                        transfer_wait_t = temp_transfer_t;
//                    if (t == 11520 && p_id == 10866) {
//                        System.out.println("transfer time: " + transfer_wait_t);
//                        System.out.println("submission time: " + requests.get(p_id).submission_t);
//
//                    }
                }
            }
            for (int p_id : route.get(k - 1).dropoff_ids) {
                if (requests.containsKey(p_id) && requests.get(p_id).transfer == true && p_id < big_N) {
                    int temp_transfer_t = requests.get(p_id).p_latest_arr_t - route.get(k - 1).exp_stop_arr_t;
                    if (temp_transfer_t > transfer_wait_t)
                        transfer_wait_t = temp_transfer_t;
                }
            }

            route.get(k).exp_stop_arr_t += transfer_wait_t;
            route.get(k).exp_stop_arr_t += dwell_time;
        }
    }

    public static int check_feasibility(Vehicle v, HashMap<Integer, Request> requests, int r_id, List<Vehicle.Route_stop> route,
                                        HashMap<Integer, List<Integer>> temp_requests, int v_cap, int base_dwell_time,
                                        int p_dwell_time, int t, int t_step, Network network, int max_wait_t) {

        int load = v.current_load;
        int feasibility = 0;

        //this part is necessary for updating temp_requests later
        // (some passengers have been picked up, so pickup time is not included in the route information)
        for (int p_id : v.active_requests) {
            temp_requests.put(p_id, Arrays.asList(requests.get(p_id).pickup_t, requests.get(p_id).dropoff_t, requests.get(p_id).submission_t));
        }

        //add the new request to the temp_requests
        if (r_id < big_N) {
            temp_requests.put(r_id, Arrays.asList(null, null, requests.get(r_id).submission_t));
        } else {
            temp_requests.put(r_id, Arrays.asList(null, null, null));
        }

        //for each stop
        for (int k = 0; k < route.size(); k++) {

            //check capacity
            load += route.get(k).pickup_ids.size() - route.get(k).dropoff_ids.size();
            if (load > v_cap) {
                feasibility = 0;
                return feasibility;
            } else {
                feasibility = 1;
            }

            //update expected arrival times
            update_exp_arr_t(k, base_dwell_time, p_dwell_time, t, t_step, network, v, route, requests, max_wait_t);

            update_temp_requests(route, temp_requests);
        }

        //check max travel time and max wait time constraints
        for (Map.Entry<Integer, List<Integer>> pair : temp_requests.entrySet()) {

                if (pair.getValue().get(0) - requests.get(pair.getKey()).submission_t > max_wait_t |
                        pair.getValue().get(1) > requests.get(pair.getKey()).p_latest_arr_t) {
                    feasibility = 0;
                    return feasibility;
                }
        }
        return feasibility;
    }

    public static void update_temp_requests(List<Vehicle.Route_stop> route, HashMap<Integer, List<Integer>> temp_requests) {

        for (int k = 0; k < route.size(); k++) {
            for (int p_id : route.get(k).pickup_ids) {
                temp_requests.put(p_id, Arrays.asList(route.get(k).exp_stop_arr_t, null, temp_requests.get(p_id).get(2)));
            }
            for (int p_id : route.get(k).dropoff_ids) {
                temp_requests.put(p_id, Arrays.asList(temp_requests.get(p_id).get(0), route.get(k).exp_stop_arr_t, temp_requests.get(p_id).get(2)));
            }
        }
    }

    public static void assignment_update(HashMap<Integer, Vehicle> vehicles, HashMap<Integer, Request> requests,
                                         Network network, int best_v_id, List<Vehicle.Route_stop> best_v_route, int r_id,
                                         Double cost, DirectedGraph graph, boolean rebalance) throws CloneNotSupportedException {

        vehicles.get(best_v_id).route_stops.clear();
        for (Vehicle.Route_stop stop : best_v_route)
            vehicles.get(best_v_id).route_stops.add(stop.clone());
//        vehicles.get(best_v_id).route_stops = best_v_route;
        if (!rebalance) {
            vehicles.get(best_v_id).active_requests.add(r_id);
        }
        vehicles.get(best_v_id).cost = cost;

        //re-find the shortest path
        List<Vehicle.Route_stop> v_stops = vehicles.get(best_v_id).route_stops;

        int next_node;
        if (vehicles.get(best_v_id).route_nodes.size() > 0) {
            next_node = vehicles.get(best_v_id).route_nodes.get(0);
            vehicles.get(best_v_id).route_nodes = new ArrayList<>();
//            vehicles.get(best_v_id).route_nodes.addAll(network.shortest_path[next_node][v_stops.get(0).stop_id]);
            vehicles.get(best_v_id).route_nodes.addAll(find_shortest_path(next_node, v_stops.get(0).stop_id, graph));
        } else {
            next_node = vehicles.get(best_v_id).current_loc;
            vehicles.get(best_v_id).route_nodes = new ArrayList<>();
//            vehicles.get(best_v_id).route_nodes.addAll(network.shortest_path[next_node][v_stops.get(0).stop_id]);
            vehicles.get(best_v_id).route_nodes.addAll(find_shortest_path(next_node, v_stops.get(0).stop_id, graph));
            //if current location and the first stop are the same, only one node is added
            //in this case, we should have the current location in the node list as well because it's a stop
//            if (vehicles.get(best_v_id).route_nodes.size()>1) {
//                vehicles.get(best_v_id).route_nodes.remove(0);
//            }
//            vehicles.get(best_v_id).tt_to_next_node = network.dist[next_node][vehicles.get(best_v_id).route_nodes.get(0)] + 30;
        }

        if (!rebalance) {
            requests.get(r_id).status = "planned";
            requests.get(r_id).assigned_veh = best_v_id;
        }

        for (int i = 0; i < v_stops.size(); i++) {
            //update route_nodes
            if (i != v_stops.size() - 1) { //since we have a pickup and a dropoff v_stop.size is always greater than 1
                //remove the last element of route_nodes since it appears again in the next shortest path
                vehicles.get(best_v_id).route_nodes.remove(vehicles.get(best_v_id).route_nodes.size() - 1);
//                vehicles.get(best_v_id).route_nodes.addAll(network.shortest_path[v_stops.get(i).stop_id][v_stops.get(i + 1).stop_id]);
                vehicles.get(best_v_id).route_nodes.addAll(find_shortest_path(v_stops.get(i).stop_id, v_stops.get(i + 1).stop_id, graph));
            }

            //update pickup times
            for (int j = 0; j < v_stops.get(i).pickup_ids.size(); j++) {
                requests.get(v_stops.get(i).pickup_ids.get(j)).pickup_t = v_stops.get(i).exp_stop_arr_t;
            }
            //update dropoff times
            for (int j = 0; j < v_stops.get(i).dropoff_ids.size(); j++) {
                requests.get(v_stops.get(i).dropoff_ids.get(j)).dropoff_t = v_stops.get(i).exp_stop_arr_t;
            }
        }
    }

    public static void print_vehicles_info(HashMap<Integer, Vehicle> vehicles, Data data) {
        for (Map.Entry<Integer, Vehicle> v : vehicles.entrySet()) {
            System.out.println("vehicle id: " + v.getKey());
            System.out.println("current load: " + v.getValue().current_load);
            System.out.println("current location: " + v.getValue().current_loc);
            System.out.println("active requests: " + v.getValue().active_requests);
            System.out.println("tt to next node: " + v.getValue().tt_to_next_node);
            System.out.println("total idle time: " + data.vehicles.get(v.getKey()).total_idle_t);
            System.out.println("route_stop: ");
            for (int i = 0; i < v.getValue().route_stops.size(); i++) {
                System.out.println("stop id: " + v.getValue().route_stops.get(i).stop_id + " pickups: " + v.getValue().route_stops.get(i).pickup_ids + " dropoffs: " + v.getValue().route_stops.get(i).dropoff_ids + " exp arr t: " + v.getValue().route_stops.get(i).exp_stop_arr_t);
            }
            System.out.println("route_node: ");
            for (int i = 0; i < v.getValue().route_nodes.size(); i++) {
                System.out.print(v.getValue().route_nodes.get(i) + ", ");
            }
            System.out.println("");
            System.out.println("cost: " + v.getValue().cost);
            System.out.println("______________________________________");
        }
    }

    public static void print_request_info(HashMap<Integer, Request> requests, int r_id) {
        System.out.println("REQUEST INFO");
        System.out.println("request id: " + r_id);
        System.out.println("origin: " + requests.get(r_id).origin);
        System.out.println("dest: " + requests.get(r_id).dest);
        System.out.println("status: " + requests.get(r_id).status);
        System.out.println("submission t: " + requests.get(r_id).submission_t);
        System.out.println("pickup t: " + requests.get(r_id).pickup_t);
        System.out.println("dropoff t: " + requests.get(r_id).dropoff_t);
        System.out.println("latest drop off time: " + requests.get(r_id).p_latest_arr_t);
        System.out.println("assigned veh: " + requests.get(r_id).assigned_veh);
        System.out.println("______________________________________");
    }

    public static void report_outputs(HashMap<Integer, Request> requests, Data data, int simulation_period,
                                      int n_vehicles, ArrayList<Float> avg_occ,
                                      int n, float theta, float beta, int n_requests,
                                      int n_stops, long start_time, CSVWriter writer, long seed, float wait_vot) {
        int sum_invehicle_time = 0;
        int sum_waiting_time = 0;
        int sum_veh_active_time = 0;
        int sum_veh_idle_time = 0;
        int sum_served_pax = 0;
        int sum_rejected = 0;
        int sum_pending = 0;
        int sum_inprocess = 0;
        int sum_submitted = 0;
        int sum_traveled_dist = 0;
        int sum_empty_traveled_dist = 0;

        for (Map.Entry<Integer, Request> r : requests.entrySet()) {
            if (r.getValue().transfer == false) {
                if (r.getValue().status.equals("alighted")) {
                    sum_invehicle_time += (r.getValue().dropoff_t - r.getValue().pickup_t);
                    sum_waiting_time += (r.getValue().pickup_t - r.getValue().submission_t);
                    sum_served_pax++;
                }
            } else {
                if (r.getValue().status.equals("alighted") && r.getKey() > big_N) {
                    sum_invehicle_time += (r.getValue().dropoff_t - requests.get(r.getKey() - big_N).pickup_t);
                    sum_waiting_time += (requests.get(r.getKey() - big_N).pickup_t - requests.get(r.getKey() - big_N).submission_t);
                    sum_served_pax++;
                }
            }
        }

        for (Map.Entry<Integer, Data.VehicleSim> v : data.vehicles.entrySet()) {
            sum_veh_active_time += (simulation_period - v.getValue().total_idle_t);
            sum_veh_idle_time += v.getValue().total_idle_t;
            sum_traveled_dist += v.getValue().total_distance;
            sum_empty_traveled_dist += v.getValue().total_empty_distance;
        }

        for (Map.Entry<Integer, Request> r : requests.entrySet()) {
            if (r.getValue().status.equals("rejected")) {
                sum_rejected += 1;
            }
        }

        for (Map.Entry<Integer, Request> r : requests.entrySet()) {
            if (r.getValue().status.equals("pending")) {
                sum_pending += 1;
            }
        }

        for (Map.Entry<Integer, Request> r : requests.entrySet()) {
            if (r.getValue().transfer == false) {
                if (r.getValue().status.equals("planned") || r.getValue().status.equals("onboard")) {
                    sum_inprocess += 1;
                }
            } else {
                if (r.getKey() > big_N && (r.getValue().status.equals("planned") || r.getValue().status.equals("onboard"))) {
                    sum_inprocess += 1;
                }
            }
        }

        float sum_occ = 0;
        for (int i = 0; i < avg_occ.size(); i++) {
            sum_occ += avg_occ.get(i);
        }

        double sum_invehicle_time_d = (double) sum_invehicle_time / 3600;
        double sum_waiting_time_d = (double) sum_waiting_time / 3600;
        double sum_veh_active_time_d = (double) sum_veh_active_time / 3600;
        double sum_veh_idle_time_d = (double) sum_veh_idle_time / 3600;
        double avg_pax_invehicle_time = (double) sum_invehicle_time / sum_served_pax / 60;
        double avg_pax_waiting_time = (double) sum_waiting_time / sum_served_pax / 60;
        double sum_traveled_dist_d = (double) sum_traveled_dist / 3600;
        double avg_traveled_dist = (double) sum_traveled_dist_d / sum_served_pax * 60;
        double sum_empty_traveled_dist_d = (double) sum_empty_traveled_dist / 3600;

        float avg_occupancy = sum_occ / (n + 1);

        DecimalFormat df = new DecimalFormat("#.##");

        System.out.println("theta: " + theta);
        System.out.println("beta: " + beta);
        System.out.println("simulation period: " + simulation_period / 3600 + " hours");
        System.out.println("number of requests: " + n_requests);
        System.out.println("number of vehicles: " + n_vehicles);
        System.out.println("total passengers' in-vehicle time (completed requests): " + df.format(sum_invehicle_time_d) + " hours");
        System.out.println("total passengers' waiting time (completed requests): " + df.format(sum_waiting_time_d) + " hours");
        System.out.println("average passengers' in-vehicle time: " + df.format(avg_pax_invehicle_time) + " minutes");
        System.out.println("average passengers' waiting time: " + df.format(avg_pax_waiting_time) + " minutes");
        System.out.println("total vehicles' active time: " + df.format(sum_veh_active_time_d) + " hours");
        System.out.println("total vehicles' idle time: " + df.format(sum_veh_idle_time_d) + " hours");
        System.out.println("average vehicles' occupancy: " + df.format(avg_occupancy));
        System.out.println("number of total passengers served: " + sum_served_pax);
        System.out.println("number of rejected requests: " + sum_rejected);
        System.out.println("number of pending requests: " + sum_pending);
        System.out.println("number of in-process requests: " + sum_inprocess);
//        System.out.println("number of submitted requests: " + sum_submitted);
        System.out.println("total traveled distance by vehicles: " + df.format(sum_traveled_dist_d) + " hours");
        System.out.println("average traveled distance by vehicles: " + df.format(avg_traveled_dist) + " minutes");
        System.out.println("total empty traveled distance by vehicles: " + df.format(sum_empty_traveled_dist_d) + " hours");

        long end_time = System.nanoTime();
        double run_time = (end_time - start_time) / 1e9;
        System.out.println("Run time: " + run_time + " sec");
        String[] line = {String.valueOf(seed), String.valueOf(simulation_period), String.valueOf(n_requests), String.valueOf(n_vehicles),
                String.valueOf(n_stops), String.valueOf(theta), String.valueOf(beta),
                String.valueOf(avg_pax_invehicle_time), String.valueOf(avg_pax_waiting_time), String.valueOf(avg_pax_invehicle_time + wait_vot * avg_pax_waiting_time),
                String.valueOf(sum_veh_idle_time_d), String.valueOf(avg_occupancy), String.valueOf(sum_served_pax),
                String.valueOf(sum_rejected + sum_pending), String.valueOf(sum_traveled_dist_d), String.valueOf(avg_traveled_dist), String.valueOf(sum_empty_traveled_dist_d), String.valueOf(run_time)};

        writer.writeNext(line);
    }

    public static void print_route(HashMap<Integer, Request> requests, HashMap<Integer, Vehicle> vehicles, int v_id, List<Vehicle.Route_stop> route) {

        System.out.println("Vehicle " + v_id);
        System.out.println("current location: " + vehicles.get(v_id).current_loc);
        System.out.println("tt to next node: " + vehicles.get(v_id).tt_to_next_node);

        if (route.size() == 0)
            route = vehicles.get(v_id).route_stops;
        List<Integer> route_nodes = vehicles.get(v_id).route_nodes;

        for (int p = 0; p < route.size(); p++) {
            System.out.println("stop id: " + route.get(p).stop_id + " exp arr t: " + route.get(p).exp_stop_arr_t +
                    " pickup ids: " + route.get(p).pickup_ids + " dropoff ids: " + route.get(p).dropoff_ids);
        }
        for (int p = 0; p < route.size(); p++) {
            for (int p_id : route.get(p).dropoff_ids) {
                if (requests.containsKey(p_id)) {
                    System.out.println("request " + p_id + " submission time: " + requests.get(p_id).submission_t);
                }
            }
        }
        System.out.println("node ids: " + route_nodes);
    }

    public static void report_requests_info(String filename, HashMap<Integer, Request> requests, HashMap<Integer, Node> nodes) throws IOException {
        //requests info output file
        File file_r = new File(filename);
        FileWriter outputfile_r = new FileWriter(file_r);
        CSVWriter writer_r = new CSVWriter(outputfile_r);
        String[] header_r = {"request_id", "submission_t", "pickup_t", "dropoff_t", "status", "pickup_x", "pickup_y", "dropoff_x", "dropoff_y"};

        writer_r.writeNext(header_r);
        for (Map.Entry<Integer, Request> r : requests.entrySet()) {
//            print_request_info(requests, r.getKey());
            String[] line = {String.valueOf(r.getKey()), String.valueOf(r.getValue().submission_t),
                    String.valueOf(r.getValue().pickup_t), String.valueOf(r.getValue().dropoff_t),
                    String.valueOf(r.getValue().status), String.valueOf(nodes.get(r.getValue().origin).x),
                    String.valueOf(nodes.get(r.getValue().origin).y), String.valueOf(nodes.get(r.getValue().dest).x),
                    String.valueOf(nodes.get(r.getValue().dest).y)};
            writer_r.writeNext(line);
        }
        writer_r.close();
    }
}