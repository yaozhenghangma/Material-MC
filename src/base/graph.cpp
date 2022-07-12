#ifndef GRAPH
#define GRAPH

#include <vector>
class Edge {
};

class Vertex{
    public:
    std::vector<Edge *> edges;
};

class DirectedEdge : public Edge {
    public:
    Vertex * vertex;
};

class UndirectedEdge : public Edge {
    public:
    Vertex * vertex1;
    Vertex * vertex2;
};

class Graph {
};

class DirectedGraph : public Graph {
    public:
    std::vector<DirectedEdge> edges;
    std::vector<Vertex> vertexes;
};

class UndirectedGraph : public Graph {
    public:
    std::vector<UndirectedEdge> edges;
    std::vector<Vertex> vertexes;
};

#endif