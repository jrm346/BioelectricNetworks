from abc import abstractmethod, ABC

from core import Cell


class Membrane(ABC):
    def __init__(self, **kargs):
        self.cell = None
        self.sites = {**kargs}

    @abstractmethod
    def bind_ligand(self, ligand: str):
        pass

    @abstractmethod
    def update(self):
        pass


class Dead(Membrane):
    """
    A Membrane for cells which are dead. It does nothing.
    """
    def __init__(self, cell: Cell):
        super().__init__()
        self.cell = cell

    def bind_ligand(self, ligand: str):
        pass

    def update(self):
        pass

