import networkx as nx
from functools import reduce


def read_data(path="/home/krzysztof/Pobrane/rosalind_grep.txt"):
    with open(path, 'r') as file:
        lines = file.read()
    return lines


def build_graph(genomes):
    graph = nx.DiGraph()
    for genome_from in genomes:
        for genome_to in genomes:
            if genome_from != genome_to:
                if all([True if cgf == cgt else False for cgf, cgt in zip(genome_from[1:], genome_to[:-1])]):
                    graph.add_edge(genome_from[1:], genome_to[1:])
    return graph


def build_circular_strings(graph, string, len, output=[]):
    node = string[1-len:]
    if node not in graph.nodes:
        return output
    for neighbor in graph[node]:
        copy_graph = nx.DiGraph(graph)
        copy_graph.remove_edge(node, neighbor)
        build_circular_strings(copy_graph, string + neighbor[-1], len, output)
        output.append(string)
    return output


data = read_data()
data = data.split('\n')
del data[-1]
print(len(data))

bruijn_graph = build_graph(data)
l = build_circular_strings(bruijn_graph, data[0], len(data[0]))
print(l)
l = [e for e in l if len(e) == len(data)-2]
print(l)
