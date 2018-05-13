from Bio import SeqIO
from graphviz import Digraph
import argparse

class Vertex:
    
    def __init__(self, seq):
        self.seq = seq
        self.coverage = 1
        self.in_edges = {}
        self.out_edges = {}
        
    def increase_coverage(self):
        self.coverage += 1

class Edge:
    
    def __init__(self,k1,k2):
        self.seq = k1 + k2[-1]
        self.coverage = 0
    
    def calc_coverage(self,c1,c2):
        self.coverage = (c1+c2)/2



class Graph:

    def __init__(self,k):
        self.vertices = {}
        self.k = k
        self.graph = Digraph()
        
    def add_read(self,read):
        read_lng = len(read)
        if read_lng < self.k:
            return
        
        kmer = read[:k]
        
        if kmer in self.vertices:
            self.vertices[kmer].increase_coverage()
        else:
            self.vertices[kmer] = Vertex(kmer)

        for next_kmer_indx in range(1,read_lng-k+1,1):
            next_kmer = read[next_kmer_indx:(next_kmer_indx+k)]
            if next_kmer in self.vertices:
                self.vertices[next_kmer].increase_coverage()
            else:
                self.vertices[next_kmer] = Vertex(next_kmer)
            
            new_edge = Edge(kmer,next_kmer)
            self.vertices[next_kmer].in_edges[kmer] = new_edge
            self.vertices[kmer].out_edges[next_kmer] = new_edge
            

            kmer = next_kmer
    
    def calc_init_edge_coverage(self):
        
        for current_vertex in self.vertices.keys():
            for next_vertex in self.vertices[current_vertex].out_edges.keys():
                self.vertices[current_vertex].out_edges[next_vertex].calc_coverage(self.vertices[current_vertex].coverage,self.vertices[next_vertex].coverage)
    
    def graphviz(self, view, outfile):

        for vertex in self.vertices:
            if view == 'full':
                self.graph.node(vertex, label=vertex)
            elif view == 'nick':
                self.graph.node(vertex, label=str(self.vertices[vertex].coverage))
            for next_vertex in self.vertices[vertex].out_edges:
                if view == 'full':
                    self.graph.edge(vertex, next_vertex, label = self.vertices[vertex].out_edges[next_vertex].seq)
                elif view == 'nick':
                    self.graph.edge(vertex, next_vertex, label = str(self.vertices[vertex].out_edges[next_vertex].coverage)+','+
                                    str(len(self.vertices[vertex].out_edges[next_vertex].seq)))
                    
        self.graph.render(filename=outfile, view=True)
        

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Graph')
    
    parser.add_argument('-f', '--file', help='fasta file', type=str, required=True)
    parser.add_argument('-k', '--kmer_size', help='kmer size', type=int, default=15)
    parser.add_argument('-v', '--view', help='choose full or nick view', type=str, required=True)
    parser.add_argument('-o', '--outfile', help='outfile with graph in .dot format', type=str, required=True)
    
    args = parser.parse_args()
    file = args.file
    k = args.kmer_size
    view = args.view
    outfile = args.outfile
    
    my_graph = Graph(k)
    
    with open(file, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            read = str(record.seq)
            my_graph.add_read(read)
            my_graph.add_read(str(record.reverse_complement().seq))

    my_graph.calc_init_edge_coverage()
    my_graph.graphviz(view, outfile)
          
