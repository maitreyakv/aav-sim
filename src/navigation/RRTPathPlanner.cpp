/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#include <lemon/list_graph.h>
#include <lemon/dijkstra.h>

#include "RRTPathPlanner.h"

// TODO: Add documentation/comments

arma::vec steer(arma::vec pos_x, arma::vec pos_y, double eta) {
    // Return the new position by moving from the first position in the direction with a small step
    return pos_x + eta * arma::normalise(pos_y - pos_x);
}

bool RRTPathPlanner::isCollide(arma::vec point_1, arma::vec point_2) {
    for (Obstacle* obstacle_ptr : this->m_obstacle_ptrs) {
        if (obstacle_ptr->isLineCollide(point_1, point_2)) { return true; }
    }
    return false;
}

std::vector<arma::vec> RRTPathPlanner::computePath(arma::vec pos_start, arma::vec pos_target) {

    // TEMP:
    double max_num_vertices = 10000;
    double eta = 0.1;
    double x_max = 3.0;
    double x_min = -3.0;
    double y_max = 3.0;
    double y_min = -3.0;
    double z_max = 3.0;
    double z_min = -3.0;

    // Get the dimension of the space
    int dim = pos_start.n_elem;

    // Get the total volume of free space
    // TEMP
    double vol_free = (x_max - x_min) * (y_max - y_min) * (z_max - z_min);

    // Get the unit sphere volume
    // TEMP
    double vol_unit_sphere = 4.0 * M_PI;

    // Compute the RRT neighborhood radius coefficient
    double gamma_star = 2.0 * pow(1.0 + 1.0 / dim, 1.0 / dim) * pow(vol_free / vol_unit_sphere, 1.0 / dim);

    // Initialize a directed graph for the RRT algorithm's tree
    lemon::ListDigraph graph;

    // Initialize a node map corresponding to the graph to assign values to the vertices in the graph
    lemon::ListDigraph::NodeMap<arma::vec> pos_map(graph);

    // Initialize an arc map corresponding to the graph to assign distances between the vertices in the graph
    lemon::ListDigraph::ArcMap<double> dist_map(graph);

    // Initialize a node map corresponding to the graph to assign costs to each vertex in the graph for RRT
    lemon::ListDigraph::NodeMap<double> cost_map(graph);

    // Add the starting node to the graph
    lemon::ListDigraph::Node start_node = graph.addNode();
    pos_map.set(start_node, pos_start);
    cost_map.set(start_node, 0.0);

    // Allocate a node for the best node near the target point and the best cost of the close points
    lemon::ListDigraph::Node best_near_target_node = start_node;
    double best_near_target_cost = std::numeric_limits<double>::infinity();

    // Perform the RRT algorithm by spawning random nodes and connecting them in the graph
    for (int n = 0; n < max_num_vertices; n++) {
        // Create a random position node in the range [0,1]
        arma::vec pos_rand = arma::zeros<arma::vec>(dim);
        for (int i = 0; i < dim; i++) { pos_rand(i) = rand() / (double) RAND_MAX; }

        // Scale the random position so it is withing the entire search domain
        // TEMP
        pos_rand(0) = (x_max - x_min) * pos_rand(0) + x_min;
        pos_rand(1) = (y_max - y_min) * pos_rand(1) + y_min;
        pos_rand(2) = (z_max - z_min) * pos_rand(2) + z_min;

        // Allocate a node for the nearest node to the sampled point and the distance
        lemon::ListDigraph::Node nearest_node;
        double nearest_distance = std::numeric_limits<double>::infinity();

        // Search through the graph to find the nearest node
        for (lemon::ListDigraph::NodeIt node(graph); node != lemon::INVALID; ++node) {
            // Compute distance between the random node and the iterated node
            double distance_to_new = arma::norm(pos_map[node] - pos_rand);

            // Replace the nearest node with the iterated node if the distance is closer
            if (distance_to_new < nearest_distance) {
                // Replace the closest node with the iterated node
                nearest_node = node;

                // Replace the closest distance
                nearest_distance = distance_to_new;
            }
        }

        // Create the new node by moving from the nearest node in the direction with step eta
        arma::vec pos_new = steer(pos_map[nearest_node], pos_rand, eta);

        // Restart the iteration if there is a collision
        if (this->isCollide(pos_map[nearest_node], pos_new)) { continue; }

        // Compute search radius for neighbors
        double radius = fmin(gamma_star * pow(log(n) / n, 1.0 / dim), eta);

        // Find the neighbors of the new node
        std::vector<lemon::ListDigraph::Node> neighbors;
        for (lemon::ListDigraph::NodeIt node(graph); node != lemon::INVALID; ++node) {
            if (arma::norm(pos_map[node] - pos_new) <= radius) { neighbors.push_back(node); }
        }

        // Add the new node to the graph
        lemon::ListDigraph::Node new_node = graph.addNode();
        pos_map.set(new_node, pos_new);

        // Set the nearest neighbor as the current best node to connect the new node to
        lemon::ListDigraph::Node best_node = nearest_node;

        // Initialize the best cost of the new node as the sum of the nearest node cost and the nearest distance
        double best_cost = cost_map[nearest_node] + nearest_distance;

        // Attempt connections along the minimum-cost path with neighbors
        for (lemon::ListDigraph::Node neighbor : neighbors) {
            if (cost_map[neighbor] + arma::norm(pos_map[neighbor] - pos_new) < best_cost
                                                                    && !this->isCollide(pos_map[neighbor], pos_new)) {
                best_node = neighbor;
                best_cost = cost_map[neighbor] + arma::norm(pos_map[neighbor] - pos_new);
            }
        }

        // Add an arc connecting the new sampled node and its best neighbor and assign the distance to the arc
        lemon::ListDigraph::Arc new_arc = graph.addArc(best_node, new_node);
        dist_map.set(new_arc, arma::norm(pos_map[best_node] - pos_new));
        cost_map.set(new_node, best_cost);

        // Attempt to rewire the tree in the neighborhood
        for (lemon::ListDigraph::Node neighbor : neighbors) {
            if (cost_map[new_node] + arma::norm(pos_new - pos_map[neighbor]) < cost_map[neighbor]
                                                                    && !this->isCollide(pos_map[neighbor], pos_new)) {
                lemon::ListDigraph::Node parent_node;
                lemon::ListDigraph::Arc parent_arc;
                parent_node = neighbor;
                for (lemon::ListDigraph::InArcIt arc(graph, neighbor); arc != lemon::INVALID; ++arc) {
                    parent_node = graph.source(arc);
                    parent_arc = arc;
                }

                if (!(parent_node == neighbor)) {
                    graph.erase(parent_arc);
                    lemon::ListDigraph::Arc new_arc = graph.addArc(new_node, neighbor);
                    dist_map.set(new_arc, arma::norm(pos_new - pos_map[neighbor]));
                    cost_map.set(neighbor, cost_map[new_node] + arma::norm(pos_new - pos_map[neighbor]));
                }
            }
        }

        // Compute distance between the target and the new node
        double distance_to_target = arma::norm(pos_map[new_node] - pos_target);

        // Replace the best node close to target with the new node if the cost is lower
        // TEMP
        if (distance_to_target < 0.1 && cost_map[new_node] < best_near_target_cost) {
            // Replace the closest node with the new node
            best_near_target_node = new_node;

            // Replace the closest distance
            best_near_target_cost = cost_map[new_node];
        }
    }

    // Allocate an array for the computed path
    std::vector<arma::vec> path;

    // Initialize a Dijkstra's algorithm based searcher to find the best path through the tree
    lemon::Dijkstra<lemon::ListDigraph, lemon::ListDigraph::ArcMap<double>, lemon::DijkstraDefaultTraits<lemon::ListDigraph,
                    lemon::ListDigraph::ArcMap<double>>> dijkstra(graph, dist_map);

    // Initialize the Dijkstra's algorithm
    dijkstra.init();

    // Use Dijkstra's algorithm to extract the shortest path from the explored tree
    dijkstra.run(start_node, best_near_target_node);

    // Starting from the target node, work backward towards the start to construct the optimal path
    for (lemon::ListDigraph::Node node = best_near_target_node; node != start_node; node = dijkstra.predNode(node)) {
        // Add the node's position to the path
        if (node != lemon::INVALID) { path.insert(path.begin(), pos_map[node]); }
    }

    // Add the initial state to the front of the path
    path.insert(path.begin(), pos_start);




    // TEMP: Testing Savitzky–Golay filtering
    for (int i = 0; i < 1000; i++) {
        for (int n = 4; n < path.size() - 4; n++) {
            path[n] = (-21.0 * path[n - 4] + 14.0 * path[n - 3] + 39.0 * path[n - 2] + 54.0 * path[n - 1]
                       + 59.0 * path[n]
                       - 21.0 * path[n + 4] + 14.0 * path[n + 3] + 39.0 * path[n + 2] + 54.0 * path[n + 1]) / 231.0;
        }
    }



    // TEMP
    for (lemon::ListDigraph::ArcIt arc(graph); arc != lemon::INVALID; ++arc) {
        arma::vec pt1 = pos_map[graph.source(arc)];
        arma::vec pt2 = pos_map[graph.target(arc)];
        std::cout << pt1(0) << "," << pt1(1) << "," << pt1(2) << "," << pt2(0) << "," << pt2(1) << "," << pt2(2) << std::endl;
    }



    // Return the computed path
    return path;
}