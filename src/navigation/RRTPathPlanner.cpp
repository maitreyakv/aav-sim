/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#include <lemon/list_graph.h>
#include <lemon/dijkstra.h>

#include "RRTPathPlanner.h"

// TODO: Add documentation/comments

std::vector<arma::vec> RRTPathPlanner::computePath(arma::vec pos_start, arma::vec pos_target) {

    // TODO: Implement me

    // TEMP:
    double max_num_vertices = 50000;
    double step = 0.1;
    double x_max = 3.0;
    double x_min = -3.0;
    double y_max = 3.0;
    double y_min = -3.0;
    double z_max = 3.0;
    double z_min = -3.0;

    // Initialize a graph for the RRT algorithm
    lemon::ListGraph graph;

    // Initialize a node map corresponding to the graph to assign values to the vertices in the graph
    lemon::ListGraph::NodeMap<arma::vec> node_map(graph);

    // Initalize and edge map corresponsing to the graph to assign distances between the vertices in the graph
    lemon::ListGraph::EdgeMap<double> edge_map(graph);

    // Add the starting node to the graph
    lemon::ListGraph::Node start_node = graph.addNode();
    node_map.set(start_node, pos_start);

    // Perform the RRT algorithm by spawning random nodes and connecting them in the graph
    for (int n = 0; n < max_num_vertices; n++) {
        // Create a random position node in the range [0,1]
        arma::vec pos_rand = {rand() / (double) RAND_MAX,
                              rand() / (double) RAND_MAX,
                              rand() / (double) RAND_MAX};

        // Scale the random position so it is withing the entire search domain
        pos_rand(0) = (x_max - x_min) * pos_rand(0) + x_min;
        pos_rand(1) = (y_max - y_min) * pos_rand(1) + y_min;
        pos_rand(2) = (z_max - z_min) * pos_rand(2) + z_min;

        // Allocate a node for the nearest node to the sampled point and the distance
        lemon::ListGraph::Node nearest_node;
        double closest_distance = std::numeric_limits<double>::infinity();

        // Search through the graph to find the nearest node
        for (lemon::ListGraph::NodeIt node(graph); node != lemon::INVALID; ++node) {
            // Compute distance between the random node and the iterated node
            double distance = arma::norm(node_map[node] - pos_rand);

            // Replace the closest node with the iterated node if the distance is closer
            if (distance < closest_distance) {
                // Replace the closest node with the iterated node
                nearest_node = node;

                // Replace the closest distance
                closest_distance = distance;
            }
        }

        // Compute the direction from the nearest node to the random node as a unit vector
        arma::vec direction = arma::normalise(pos_rand - node_map[nearest_node]);

        // Create the new node by moving from the nearest node in the direction with a small step
        arma::vec pos_new = node_map[nearest_node] + step * direction;

        // Create the new node in the graph
        lemon::ListGraph::Node new_node = graph.addNode();
        node_map.set(new_node, pos_new);

        // Add an edge connecting the new sampled node and its nearest neighbor and assign the distance to the edge
        lemon::ListGraph::Edge new_edge = graph.addEdge(new_node, nearest_node);
        edge_map.set(new_edge, step);
    }

    // Allocate a node for the nearest node to the target point and the distance
    lemon::ListGraph::Node nearest_node;
    double closest_distance = std::numeric_limits<double>::infinity();

    // Search through the graph to find the nearest node to the target
    for (lemon::ListGraph::NodeIt node(graph); node != lemon::INVALID; ++node) {
        // Compute distance between the target and the iterated node
        double distance = arma::norm(node_map[node] - pos_target);

        // Replace the closest node with the iterated node if the distance is closer
        if (distance < closest_distance) {
            // Replace the closest node with the iterated node
            nearest_node = node;

            // Replace the closest distance
            closest_distance = distance;
        }
    }

    // Allocate an array for the computed path
    std::vector<arma::vec> path;

    // Initialize a Dijkstra's algorithm based searcher to find the best path through the tree
    lemon::Dijkstra<lemon::ListGraph, lemon::ListGraph::EdgeMap<double>, lemon::DijkstraDefaultTraits<lemon::ListGraph,
                    lemon::ListGraph::EdgeMap<double>>> dijkstra(graph, edge_map);

    // Initialize the Dijkstra's algorithm
    dijkstra.init();

    // Use Dijkstra's algorithm to extract the shortest path from the explored tree
    dijkstra.run(start_node, nearest_node);

    // Starting from the target node, work backward towards the start to construct the optimal path
    for (lemon::ListGraph::Node node = nearest_node; node != start_node; node = dijkstra.predNode(node)) {
        // Add the node's position to the path
        if (node != lemon::INVALID) { path.insert(path.begin(), node_map[node]); }
    }

    // Add the initial state to the front of the path
    path.insert(path.begin(), pos_start);

    // TEMP
    /**
    for (lemon::ListGraph::EdgeIt edge(graph); edge != lemon::INVALID; ++edge) {
        arma::vec pt1 = node_map[graph.u(edge)];
        arma::vec pt2 = node_map[graph.v(edge)];
        std::cout << pt1(0) << "," << pt1(1) << "," << pt1(2) << "," << pt2(0) << "," << pt2(1) << "," << pt2(2) << std::endl;
    }
     */

    // Return the computed path
    return path;
}