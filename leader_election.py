import csv

import pandas as pd
from random import choice

from sklearn.linear_model import LinearRegression
from math import floor, log

from benetwork import BENetwork
from bioelectric_events import BioelectricEvent, SimpleLigandEvent
from becell import BECell
from ligand_sites import SimpleLigandSite, LigandSite


class WinnerEvent(BioelectricEvent):
    """
    The bioelectric event is triggered when a cell wins leader election. It triggers a morphic change in the
    cell and sends a ligand to the adjacent cells that causes a morphic change in them to looser cells if they have the
    `LeaderElection` membrane
    """

    def __init__(self):
        super().__init__(lambda x: True if x >= 2 else False)

    def update(self, cell: BECell):
        if self.firing_function(cell.potential):
            for node in cell.edges:
                node.bind_ligand("vic")
            cell.morphic_tansform(0, 0, 0, [], [], "winner")


class CompetingEvent(SimpleLigandEvent):
    """
    The bioelectric event for cells still in competition. If triggered it sends a ligand to all adjacent nodes and has
    a positive offset of 0.5
    """
    def __init__(self):
        def ff(potential):
            if potential < 0.5:
                return False
            elif potential >= 1:
                return True
            else:
                return choice((True, False))
        super().__init__(ff, "adj", 0.5)


class WinnerSite(LigandSite):
    def __init__(self):

        super().__init__("vic", 1)

    def update(self, cell: BECell):
        if self.ligand_count > 0:
            cell.morphic_tansform(0, 0, 0, [], [], "looser")
            self.ligand_count = 0


class CompetitionSite(SimpleLigandSite):
    def __init__(self, sites: int):
        super().__init__("adj", sites, -3/2)


class LeaderElection(BECell):
    def __init__(self, maximum_degree: int, competition_sites: int):
        super().__init__(0, -2, 2, 0.5, {CompetingEvent(), WinnerEvent()},
                         {CompetitionSite(competition_sites), WinnerSite()},
                         maximum_degree, "Competing")


class LeaderElectionNetwork(BENetwork):

    def __init__(self, number_of_cells: int, max_degree: int, number_of_edges: int, competition_sites: int = 1):
        super().__init__(number_of_cells, max_degree, number_of_edges)
        self.competition_sites = competition_sites
        built = False
        while not built:
            try:
                self.build_network()
                built = True
            except Exception:
                pass

    def build_network(self):
        self.cells = tuple(LeaderElection(self.max_degree, self.competition_sites) for _ in range(self.number_of_cells))
        connected = set()
        unconnected = list(self.cells)
        for _ in range(self.number_of_edges + 1):
            if not connected:
                connected.add(unconnected.pop())
            elif unconnected:
                cell1 = unconnected.pop()
                cell2 = choice(tuple(connected))
                cell2.add_edge(cell1)
                connected.add(cell1)
                if cell2.maxed_out():
                    connected.remove(cell2)
            else:
                cell1 = choice(tuple(connected))
                cell2 = choice(tuple(connected - cell1.edges - {cell1}))
                cell2.add_edge(cell1)
                if cell2.maxed_out():
                    connected.remove(cell2)
                if cell1.maxed_out():
                    connected.remove(cell1)

    def stopping_condition(self) -> bool:
        return all(cell.status in {"looser", "winner"} for cell in self.cells)


def score(tup):
    i, j, k, l = tup
    scores = []
    for _ in range(1000):
        net = LeaderElectionNetwork(i, j, k, l)
        scores.append(net.run_network())
    return i, j, k, max(scores), sum(scores)/1000


if __name__ == "__main__":
    results = []
    configs = []
    for i in range(1, 11):
        for j in range(2, i):
            for k in range(i-1, floor(i*j/2)):
                configs.append((i, j, k, 1))

    from multiprocessing import Pool
    p = Pool(4)
    results = p.map(score, configs)

    with open("results.csv", "w") as outfile:
        writer = csv.writer(outfile)
        outfile.writelines(results)
