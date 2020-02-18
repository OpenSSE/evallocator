use std::collections::VecDeque;

#[derive(Debug, Clone)]
struct Edge {
    pub label: u64,
    pub start: usize, // pointer to the edge's start
    pub end: usize,   // pointer to the edge's end
    pub capacity: u64,
}

#[derive(Debug)]
struct Vertex {
    pub label: u64,
    pub in_edges: Vec<usize>,
    pub out_edges: Vec<usize>,
}

#[derive(Debug)]
struct Graph {
    vertices: Vec<Vertex>,
    edges: Vec<Edge>,
}

#[derive(Debug)]
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

    // fn add_vertex(&mut self, v: Vertex) -> usize {
    //     self.vertices.push(v);
    //     self.vertices.len() - 1
    // }

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
                capacity: cap,
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
    // source to the sink
    fn find_path(&self, source: usize, sink: usize, alg: TraversalAlgorithm) -> Option<Vec<usize>> {
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

        if visited[sink] {
            // it is faster to compute the size of the vector, and then
            // construct it, than it is to rely on multiple dynamic allocations

            let mut v = sink;
            let mut size = 0;
            loop {
                match parent_edge[v] {
                    None => break,
                    Some(e) => {
                        v = self.edges[e].start;
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

            Some(path)
        } else {
            None
        }
    }

    fn find_path_bfs(&self, source: usize, sink: usize) -> Option<Vec<usize>> {
        self.find_path(source, sink, TraversalAlgorithm::BreadthFirstSearch)
    }

    fn find_path_dfs(&self, source: usize, sink: usize) -> Option<Vec<usize>> {
        self.find_path(source, sink, TraversalAlgorithm::DepthFirstSearch)
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Max Flow");

    let mut g = Graph::new();

    let a = g.add_vertex(31);
    let b = g.add_vertex(32);
    let c = g.add_vertex(33);
    let d = g.add_vertex(34);
    let e = g.add_vertex(35);
    let f = g.add_vertex(36);

    g.add_edge(40, a, b, 10);
    g.add_edge(41, b, c, 5);
    g.add_edge(43, b, d, 7);
    g.add_edge(44, c, e, 7);
    g.add_edge(45, e, f, 7);
    g.add_edge(46, d, f, 7);

    println!("Graph: {:?}", g);

    // println!("Path a->c : {:?}", g.find_path(a, c));
    // println!("Path a->d : {:?}", g.find_path(a, d));
    println!("Path a->f : {:?}", g.find_path_bfs(a, f));
    println!("Path a->f : {:?}", g.find_path_dfs(a, f));

    Ok(())
}
