use std::collections::VecDeque;
extern crate rand;
use rand::prelude::*;

#[derive(Debug, Clone)]
struct Edge {
    pub label: u64,
    pub start: usize, // pointer to the edge's start
    pub end: usize,   // pointer to the edge's end
    pub capacity: i64,
}

#[derive(Debug, Clone)]
struct Vertex {
    pub label: u64,
    pub in_edges: Vec<usize>,
    pub out_edges: Vec<usize>,
}

#[derive(Debug, Clone)]
struct Graph {
    vertices: Vec<Vertex>,
    edges: Vec<Edge>,
}

#[derive(Debug, Clone, Copy)]
enum TraversalAlgorithm {
    DepthFirstSearch,
    BreadthFirstSearch,
}

impl Graph {
    fn new() -> Graph {
        Graph {
            vertices: Vec::new(),
            edges: Vec::new(),
        }
    }

    fn new_with_vertices(n_vertices: usize) -> Graph {
        let mut vertices: Vec<Vertex> = Vec::with_capacity(n_vertices);
        for i in 0..n_vertices {
            let v = Vertex {
                label: i as u64,
                in_edges: Vec::new(),
                out_edges: Vec::new(),
            };
            vertices.push(v);
        }
        Graph {
            vertices: vertices,
            edges: Vec::new(),
        }
    }

    // Construct a residual graph for the Ford-Fulkerson algorithm
    // from the input graph
    // An important invariant is that an edge with index e in the original
    // graph correspond the the edge with index e in the residual graph and,
    // the reversed edge in the residual graph has index 2e.
    fn new_residual_graph(graph: &Graph) -> Graph {
        let mut res_graph = graph.clone();
        let mut inv_edges = res_graph.edges.clone();

        let n_edges = inv_edges.len();

        // create inverted edges
        for e in inv_edges.iter_mut() {
            let v = e.end;
            e.end = e.start;
            e.start = v;
            e.capacity = 0; // these newly created edges have no real capacity
        }

        // and append them to the regular edges
        res_graph.edges.append(&mut inv_edges);

        for v in res_graph.vertices.iter_mut() {
            let mut reversed_out: Vec<usize> = v.out_edges.iter().map(|e| e + n_edges).collect();
            let mut reversed_in: Vec<usize> = v.in_edges.iter().map(|e| e + n_edges).collect();
            v.out_edges.append(&mut reversed_in);
            v.in_edges.append(&mut reversed_out);
        }

        res_graph
    }

    fn add_vertex(&mut self, l: u64) -> usize {
        self.vertices.push(Vertex {
            label: l,
            in_edges: vec![],
            out_edges: vec![],
        });
        self.vertices.len() - 1
    }

    fn add_edge(&mut self, label: u64, start: usize, end: usize, cap: u64) {
        if start < self.vertices.len() && end < self.vertices.len() {
            let e = Edge {
                label,
                start,
                end,
                capacity: cap as i64,
            };
            self.edges.push(e);
            let e_index = self.edges.len() - 1;
            self.vertices[start].out_edges.push(e_index);
            self.vertices[end].in_edges.push(e_index);
        } else if start >= self.vertices.len() {
            panic!("Invalid starting vertex");
        } else {
            panic!("Invalid ending vertex");
        }
    }

    // Takes a source and a sink and finds a path from the source to the sink.
    // The path is (optionally) returned as the list of edges leading from the
    // source to the sink, together with the capacity of the path
    fn find_path(
        &self,
        source: usize,
        sink: usize,
        alg: TraversalAlgorithm,
    ) -> Option<(Vec<usize>, i64)> {
        if source >= self.vertices.len() {
            panic!("Invalid source");
        }
        if sink >= self.vertices.len() {
            panic!("Invalid sink");
        }
        let mut visited = vec![false; self.vertices.len()];
        let mut parent_edge: Vec<Option<usize>> = vec![None; self.vertices.len()];
        let mut queue = VecDeque::<usize>::with_capacity(self.vertices.len());

        queue.push_front(source);
        visited[source] = true;

        // println!("Source {}", source);
        // println!("Sink {}", sink);

        let mut found_sink = false; // flag for early exit
                                    // as we have two nested loops (and we do not want to use gotos -- this
                                    // is Rust after all), this flag is necessary
        while !found_sink {
            match queue.pop_front() {
                None => {
                    // println!("\nEmpty queue");
                    break;
                }
                Some(v) => {
                    // println!("\nVertex {}", v);

                    let out_edges_index = &(self.vertices[v].out_edges);

                    for e in out_edges_index {
                        // println!("Edge {}", *e);
                        if self.edges[*e].capacity > 0 {
                            let dest = self.edges[*e].end;
                            if !visited[dest] {
                                match alg {
                                    TraversalAlgorithm::BreadthFirstSearch => queue.push_back(dest),
                                    TraversalAlgorithm::DepthFirstSearch => queue.push_front(dest),
                                }
                                parent_edge[dest] = Some(*e);
                                visited[dest] = true;
                                if dest == sink {
                                    found_sink = true;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

        if visited[sink] {
            // it is faster to compute the size of the vector, and then
            // construct it, than it is to rely on multiple dynamic allocations
            let mut capacity = std::i64::MAX;
            let mut v = sink;
            let mut size = 0;
            loop {
                match parent_edge[v] {
                    None => break,
                    Some(e) => {
                        v = self.edges[e].start;
                        capacity = capacity.min(self.edges[e].capacity);
                        size += 1;
                    }
                }
            }
            assert_ne!(size, 0);
            let mut path = vec![0; size];

            v = sink;
            let mut i = 0;

            loop {
                match parent_edge[v] {
                    None => break,
                    Some(e) => {
                        path[size - 1 - i] = e;
                        v = self.edges[e].start;

                        i += 1;
                    }
                }
            }

            Some((path, capacity))
        } else {
            None
        }
    }

    fn find_path_bfs(&self, source: usize, sink: usize) -> Option<(Vec<usize>, i64)> {
        self.find_path(source, sink, TraversalAlgorithm::BreadthFirstSearch)
    }

    fn find_path_dfs(&self, source: usize, sink: usize) -> Option<(Vec<usize>, i64)> {
        self.find_path(source, sink, TraversalAlgorithm::DepthFirstSearch)
    }

    fn compute_max_flow(
        &self,
        source: usize,
        sink: usize,
        traversal_alg: TraversalAlgorithm,
    ) -> Graph {
        let mut res_graph = Graph::new_residual_graph(self);

        let mut max_flow = 0;
        let n_original_edges = self.edges.len();
        loop {
            match res_graph.find_path(source, sink, traversal_alg) {
                None => break,
                Some((path, path_flow)) => {
                    for e in path {
                        // update capacities
                        res_graph.edges[e].capacity -= path_flow;
                        // get the reversed edge
                        let rev_edge = if e < n_original_edges {
                            e + n_original_edges
                        } else {
                            e - n_original_edges
                        };
                        res_graph.edges[rev_edge].capacity += path_flow;
                    }
                    max_flow += path_flow;
                }
            }
        }

        println!("Max flow {}", max_flow);

        res_graph.transform_residual_to_flow();
        res_graph
    }

    fn transform_residual_to_flow(&mut self) {
        assert_eq!(self.edges.len() % 2, 0);

        let n_real_edges = self.edges.len() / 2;

        for v in self.vertices.iter_mut() {
            v.in_edges = v
                .in_edges
                .iter()
                .filter(|e| n_real_edges > **e)
                .map(|e| *e) // why do we need that ??
                .collect();

            v.out_edges = v
                .out_edges
                .iter()
                .filter(|e| n_real_edges > **e)
                .map(|e| *e)
                .collect();
        }

        // we would like to do this, but due to double mutable borrow in the closure, we have an issue
        // self.edges[..n_real_edges]
        //     .iter_mut()
        //     .enumerate()
        //     .for_each(|(i, mut e)| {
        //         e.capacity = self.edges[i + n_real_edges].capacity;
        //     });

        for i in 0..n_real_edges {
            self.edges[i].capacity = self.edges[i + n_real_edges].capacity;
        }

        // remove the now unnecessary edges
        self.edges.truncate(n_real_edges);
    }
}

fn generate_random_graph(n_vertices: usize, n_edges: usize) -> Graph {
    let mut graph = Graph::new();

    for i in 0..n_vertices {
        graph.add_vertex(i as u64);
    }

    let mut rng = thread_rng();

    for i in 0..n_edges {
        let s = rng.gen_range(0, n_vertices);
        let mut e = rng.gen_range(0, n_vertices);

        while e == s {
            // avoid loops
            e = rng.gen_range(0, n_vertices);
        }
        graph.add_edge(i as u64, s, e, 1);
    }

    graph
}

fn flow_alloc(n: usize, m: usize, max_len: usize) -> Vec<usize> {
    let mut remaining_elements = n;

    let mut rng = thread_rng();

    // create a new graph with m+2 vertices: one per bucket + a source and a sink
    // by convention, the source has index m and the sink has index m+1
    let mut graph = Graph::new_with_vertices(m + 2);

    let list_index: u64 = 0;

    while remaining_elements != 0 {
        let l: usize = rng.gen_range(0, max_len.min(remaining_elements)) + 1;
        let h1: usize = rng.gen_range(0, m);
        let h2: usize = rng.gen_range(0, m);

        // Adding a list of size l consists in adding an edge of capacity l
        // between to random vertices
        g.add_vertex(list_index, h1, h2, l);

        remaining_elements -= l;
        list_index += 1;
    }

    // OK, now we have to add the edges originating from the source and the
    // ones ending at the sink to 'encode' overflowing and underflowing nodes.

    // It is time to max flow !
    let ff = g.compute_max_flow(0, n_vertices - 1, TraversalAlgorithm::DepthFirstSearch);

    // and now, look at the results.
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Max Flow");

    // let mut g = Graph::new();

    // let a = g.add_vertex(31);
    // let b = g.add_vertex(32);
    // let c = g.add_vertex(33);
    // let d = g.add_vertex(34);
    // let e = g.add_vertex(35);
    // let f = g.add_vertex(36);

    // g.add_edge(40, a, b, 10);
    // g.add_edge(41, b, c, 5);
    // g.add_edge(43, b, d, 7);
    // g.add_edge(44, c, e, 7);
    // g.add_edge(45, e, f, 7);
    // g.add_edge(46, d, f, 7);

    let n_edges = 1 << 20 as usize;
    let n_vertices = 1 << 10 as usize;

    println!(
        "Generate graph with {} vertices and {} edges of capacity 1",
        n_edges, n_vertices
    );

    let g = generate_random_graph(n_vertices, n_edges);
    // println!("Graph: {:?}", g);
    println!("Graph generated!");

    // println!("Path a->c : {:?}", g.find_path(a, c));
    // println!("Path a->d : {:?}", g.find_path(a, d));
    // println!("Path a->f : {:?}", g.find_path_bfs(a, f));
    // println!("Path a->f : {:?}", g.find_path_dfs(a, f));

    println!("Start computing max flow...");
    let ff = g.compute_max_flow(0, n_vertices - 1, TraversalAlgorithm::DepthFirstSearch);
    println!("Max flow graph computed...");
    // println!("\n\nFF Graph: {:?}", ff);

    Ok(())
}
