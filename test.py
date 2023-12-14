class Node:
    """ Node in a de Bruijn graph, representing a k-1 mer. """

    def __init__(self, km1mer):
        self.km1mer = km1mer
        self.nin = 0
        self.nout = 0

    def is_semi_balanced(self):
        return abs(self.nin - self.nout) == 1

    def is_balanced(self):
        return self.nin == self.nout

    # def __hash__(self):
    #     return hash(str(self.km1mer))

    def __str__(self):
        return str(self.km1mer)


class DeBruijnGraph:
    """ A de Bruijn multigraph built from a collection of strings. """

    def __init__(self, str_iter, k):
        """ Build de Bruijn multigraph given string iterator and k-mer length k """
        self.G = {}     # multimap from nodes to neighbors
        self.nodes = {}  # maps k-1-mers to Node objects
        self.head = None
        self.tail = None
        self.nsemi = 0
        self.nbal = 0
        self.nneither = 0

        self._build_graph(str_iter, k)

    def _build_graph(self, str_iter, k):
        for st in str_iter:
            for kmer, km1L, km1R in self.chop(st, k):
                nodeL = self._get_or_create_node(km1L)
                nodeR = self._get_or_create_node(km1R)
                nodeL.nout += 1
                nodeR.nin += 1
                self.G.setdefault(nodeL, []).append(nodeR)

        for node in self.nodes.values():
            if node.is_balanced():
                self.nbal += 1
            elif node.is_semi_balanced():
                if node.nin == node.nout + 1:
                    self.tail = node
                if node.nin == node.nout - 1:
                    self.head = node
                self.nsemi += 1
            else:
                self.nneither += 1

    def _get_or_create_node(self, km1mer):
        if km1mer in self.nodes:
            return self.nodes[km1mer]
        else:
            new_node = Node(km1mer)
            self.nodes[km1mer] = new_node
            return new_node

    def chop(self, st, k):
        """ Chop a string up into k mers of given length """
        result = []
        for i in range(len(st)-(k-1)):
            result.append((st[i:i+k], st[i:i+k-1], st[i+1:i+k]))
        return result

    def nnodes(self):
        """ Return # nodes """
        return len(self.nodes)

    def nedges(self):
        """ Return # edges """
        return len(self.G)

    def has_eulerian_path(self):
        """ Return true iff graph has Eulerian path. """
        return self.nneither == 0 and self.nsemi == 2

    def has_eulerian_cycle(self):
        """ Return true iff graph has Eulerian cycle. """
        return self.nneither == 0 and self.nsemi == 0

    def is_eulerian(self):
        """ Return true iff graph has Eulerian path or cycle """
        return self.has_eulerian_path() or self.has_eulerian_cycle()

    def eulerian_path(self):
        """ Find and return Eulerian path or cycle (as appropriate) """
        assert self.is_eulerian()
        g = self.G

        if self.has_eulerian_path():
            g = g.copy()
            assert self.head is not None
            assert self.tail is not None
            g.setdefault(self.tail, []).append(self.head)

        tour = []
        src = list(g.keys())[0]  # pick arbitrary starting node

        def _visit(n):
            while len(g[n]) > 0:
                dst = g[n].pop()
                _visit(dst)
            tour.append(n)

        _visit(src)
        tour = tour[::-1][:-1]

        if self.has_eulerian_path():
            sti = tour.index(self.head)
            tour = tour[sti:] + tour[:sti]

        return list(map(str, tour))

    def assemble_genome(self, tour):
        assert self.is_eulerian()
        s = ''
        for i, seq in enumerate(tour):
            if i == 0:
                s += seq
            else:
                s += seq[-1]
        return s

# sequences = ["GTAA", "AAGC", "TCGT"]
sequences = ['CTCTAGCC', 'TAGCCCCCT']


k = 6

graph = DeBruijnGraph(sequences, k)

for key, values in graph.G.items():
    print(key, list(map(str, values)))

# for key, values in graph.nodes.items():
#     print(key, values)

print(graph.G)

tour = graph.eulerian_path()

print(graph.G)

genome = graph.assemble_genome(tour)

print(genome)