import networkx as nx


def read_data(path="/home/krzysztof/Pobrane/rosalind_grep.txt"):
    with open(path, 'r') as file:
        lines = file.read()
    return lines


def build_graph(genomes):
    graph = nx.MultiDiGraph()
    for genome in genomes[1:]:
        graph.add_edge(genome[:-1], genome[1:])
    return graph


def build_circular_strings(graph, string, len, output=[]):
    node = string[1-len:]
    for neighbor in graph[node]:
        copy_graph = nx.MultiDiGraph(graph)
        copy_graph.remove_edge(node, neighbor)
        build_circular_strings(copy_graph, string + neighbor[-1], len, output)
        output.append(string)
    return output


data = read_data()
data = data.split('\n')

bruijn_graph = build_graph(data)
cycles = build_circular_strings(bruijn_graph, data[0], len(data[0]))
cycles = [cycle for cycle in cycles if len(cycle) == len(data)]
for cycle in cycles:
    print(cycle[:-1])
