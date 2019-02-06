import networkx as nx


def read_data(path="/home/krzysztof/Pobrane/rosalind_grep.txt"):
    with open(path, 'r') as file:
        lines = file.read()
    return lines


def build_graph(genomes):
    """
    Buduje graf de Bruijna. W wierzchołkach grafu znajdują się ciągi znakowe o jeden mniejsze niż fragmenty DNA.
    :param genomes: Tablica fragmentów genotypu.
    :return: Graf de Bruijna.
    """
    graph = nx.MultiDiGraph()
    for genome in genomes[1:]:
        graph.add_edge(genome[:-1], genome[1:])
    return graph


def build_circular_strings(graph, string, len, output=[]):
    """
    Znajduje wszystkie unikalne sekwencje, które można złożyć z podanych fragmentów. W tym celu odnajduje
    wszystkie cykle Eulera występujące w grafie. Jako że liczba krawędzi między wierzchołkami jest nie wielka,
    do tego celu wykorzystuje metodę brote force.
    :param graph: Graf de Bruijna, przechowywany w klasie nx.MultiDiGraph
    :param string: Inicjijemy fragmentem sekwencji, który powinien znaleść się na początku całej sekwencji.
    :param len: Długość fragmentów genotypu.
    :param output: Nie powinno się podawać argumentu. Służy do rekurencji.
    :return: Tablica wszystkich zwierająca sekwencje uzyskane, przy sprawdzaniu wszystkich cykli w grafie. Żeby uzykać
    pełne ciągi należy wybrać tylko te równe liczbie podanych sekwencji.
    """
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
