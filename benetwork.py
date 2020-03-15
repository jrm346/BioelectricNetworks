from __future__ import annotations
from abc import ABC, abstractmethod

import math


class BENetwork(ABC):
    @abstractmethod
    def __init__(self, number_of_cells: int, max_degree: int, number_of_edges: int):
        self.number_of_cells = number_of_cells
        self.max_degree = max_degree
        self.number_of_edges = number_of_edges
        if self.number_of_edges > math.floor((self.number_of_cells * self.max_degree)/2):
            raise ValueError(f"Maximum number of edges for a graph of {self.number_of_cells} cells with maximum degree "
                             f"{self.max_degree} is {math.floor((self.number_of_cells * self.max_degree)/2)}")
        self.cells = []
        self.rounds = 0

    def simulate_round(self):
        for cell in self.cells:
            cell.bioelectric_updates()

        for cell in self.cells:
            cell.ligand_updates()
            cell.gradient_update()
            cell.check_min_potential()
            cell.update_potential()

        self.rounds += 1

    @abstractmethod
    def build_network(self):
        pass

    @abstractmethod
    def stopping_condition(self) -> bool:
        pass

    def run_network(self) -> int:
        while not self.stopping_condition():
            self.simulate_round()
        return self.rounds
