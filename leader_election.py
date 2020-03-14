from collections import Counter
from random import choice

from bioelectric_events import BioelectricEvent, SimpleLigandEvent
from core import Cell
from membranes import Dead, Membrane


class WinnerEvent(BioelectricEvent):
    """
    The bioelectric event is triggered when a cell wins leader election. It triggers a morphic change in the
    cell and sends a ligand to the adjacent cells that causes a morphic change in them to looser cells if they have the
    `LeaderElection` membrane
    """

    def __init__(self):
        super().__init__(lambda x: True if x >= 2 else False)

    def __call__(self, cell: Cell):
        if self.firing_function(cell.potential):
            for node in cell.edges:
                node.bind_ligand("vic")
            cell.morphic_tansform(0, 0, 0, [], Dead, "winner")


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
                choice((True, False))
        super().__init__(ff, "adj", 0.5)


class LeaderMembrane(Membrane):
    """
    The membrane for cells competing in leader election

    :param adj: How many sites available for ligands transmitted by cells still in competition.
    """
    def __init__(self, adj=1):
        super().__init__(**locals(), vic=1)
        self.ligands = Counter()

    def bind_ligand(self, ligand: str):
        """
        Bind a ligand if there is an available site

        :param ligand: the ligand to bind
        :return: None
        """
        if self.ligands[ligand] < self.sites[ligand]:
            self.ligands[ligand] += 1
        else:
            pass

    def update(self):
        """
        Updates the cell the membrane belongs to. If a neighbor fired decreases potential by -3/2. If a ligand from a
        winner arrives triggers a morphic change to a looser

        :return: None
        """
        if self.ligands["vic"] > 0:
            self.cell.morphic_tansform(0, 0, 0, [], Dead(self.cell), "looser")
        elif self.ligands["adj"] > 0:
            self.cell.potential_change += -3/2
        self.ligands.clear()


class LeaderElection(Cell):
    def __init__(self, maximum_degree: int, ligand_sites={"adj": 1}):
        super().__init__(0, -2, 2, 0.5, {CompetingEvent(), WinnerEvent()}, LeaderMembrane, maximum_degree,
                         "Competing", ligand_sites)
