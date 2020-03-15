from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Callable

from becell import BECell


class BioelectricEvent(ABC):
    def __init__(self, firing_function: Callable[[float], bool], *args):
        self.firing_function = firing_function

    @abstractmethod
    def update(self, cell: BECell):
        """
        The bioelectric event. Takes a potential and checks if the event fires.

        :param cell: the cell the event is to be applied to
        :return: None
        """
        pass


class SimpleLigandEvent(BioelectricEvent):
    def __init__(self, firing_function: Callable[[float], bool], ligand: str, offset: float):
        super().__init__(firing_function)
        self.ligand = ligand
        self.offset = offset

    def update(self, cell: BECell):
        """
        The bioelectric event. Takes a potential and checks if the event fires. If it does

        :param cell: the cell the event is to be applied to
        :return: None
        """
        if self.firing_function(cell.potential):
            for node in cell.edges:
                node.bind_ligand(self.ligand)

            cell.potential_change += self.offset


