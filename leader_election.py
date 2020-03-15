from random import choice

from bioelectric_events import BioelectricEvent, SimpleLigandEvent
from core import Cell
from ligand_sites import SimpleLigandSite, LigandSite


class WinnerEvent(BioelectricEvent):
    """
    The bioelectric event is triggered when a cell wins leader election. It triggers a morphic change in the
    cell and sends a ligand to the adjacent cells that causes a morphic change in them to looser cells if they have the
    `LeaderElection` membrane
    """

    def __init__(self):
        super().__init__(lambda x: True if x >= 2 else False)

    def update(self, cell: Cell):
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
                choice((True, False))
        super().__init__(ff, "adj", 0.5)


class WinnerSite(LigandSite):
    def __init__(self):

        super().__init__("vic", 1)

    def update(self, cell: Cell):
        if self.ligand_count > 0:
            cell.morphic_tansform(0, 0, 0, [], [], "looser")
            self.ligand_count = 0


class CompetitionSite(SimpleLigandSite):
    def __init__(self, sites: int):
        super().__init__("adj", sites, -3/2)


class LeaderElection(Cell):
    def __init__(self, maximum_degree: int, competition_sites: int):
        super().__init__(0, -2, 2, 0.5, {CompetingEvent(), WinnerEvent()},
                         {CompetitionSite(competition_sites), WinnerSite()},
                         maximum_degree, "Competing")
