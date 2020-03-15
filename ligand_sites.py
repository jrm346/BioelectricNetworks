from abc import abstractmethod, ABC

from core import Cell


class LigandSite:
    def __init__(self, ligand: str, sites: int):
        self.ligand = ligand
        self.sites = sites
        self.ligand_count = 0

    @abstractmethod
    def update(self, cell: Cell):
        pass

    def bind(self):
        if self.ligand_count < self.sites:
            self.ligand_count += 1


class SimpleLigandSite(LigandSite):
    def __init__(self, ligand: str, sites: int, potential_update: float):
        super().__init__(ligand, sites)
        self.potential_update = potential_update

    def update(self, cell: Cell):
        if self.ligand_count > 0:
            cell.potential_change += self.potential_update * self.ligand_count
            self.ligand_count = 0
