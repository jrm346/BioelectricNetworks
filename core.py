from __future__ import annotations

from typing import Iterable
from bioelectric_events import BioelectricEvent
from ligand_sites import LigandSite


class Cell:
    def __init__(self, initial_potential: float, minimum_potential: float,
                 equilibrium_potential: float, equilibrium_gradient: float,
                 bioelectric_events: Iterable[BioelectricEvent], ligand_sites: Iterable[LigandSite],
                 maximum_degree: int, status: str):
        self.type = status
        self.potential = initial_potential
        self.equilibrium_potential = equilibrium_potential
        self.equilibrium_gradient = equilibrium_gradient
        self.minimum_potential = minimum_potential
        self.bioelectric_events = bioelectric_events
        self.maximum_degree = maximum_degree
        self.edges = set()
        self.ligand_sites = {site.ligand: site for site in ligand_sites}
        self.potential_change = 0

    def add_edge(self, node: Cell) -> None:
        """
        Add an edge to the given node. Edges are undirected so both nodes are updated in place.

        :param node: The node you want to add an edge to.
        :return: None
        """
        if len(self.edges) >= self.maximum_degree:
            raise ValueError(f"Node: {self} already has the maximum number of edges")
        if len(node.edges) >= node.maximum_degree:
            raise ValueError(f"Node: {node} already has the maximum number of edges")
        if node in self.edges:
            raise ValueError(f"Node: {self} already has an edge to node: {node}")
        self.edges.add(node)
        node.edges.add(self)

    def bind_ligand(self, ligand: str) -> None:
        """
        A function for receiving ligands.

        :param ligand: The ligand that is sent
        :return: None
        """
        if ligand in self.ligand_sites:
            self.ligand_sites[ligand].bind()

    def attempt_events(self) -> None:
        """
        Attempts to perform the bioelectric event.

        :return: None
        """
        for event in self.bioelectric_events:
            event.update(self)

    def apply_membrane_function(self) -> None:
        """
        Apply the membrane function to the received ligands, (often clears the ligands and updates the
        potential_change).

        :return: None
        """
        for site in self.ligand_sites.values():
            site.update(self)

    def apply_gradient(self) -> None:
        """
        Apply the update based on equilibrium gradient

        :return: None
        """
        z = self.potential - self.equilibrium_potential
        change = 0
        if z >= self.equilibrium_gradient:
            change = -self.equilibrium_gradient
        elif 0 < z < self.equilibrium_gradient:
            change = -z
        elif z == 0:
            change = 0
        elif -self.equilibrium_gradient < z < 0:
            change = z
        elif z < -self.equilibrium_gradient:
            change = self.equilibrium_gradient

        self.potential_change += change

    def check_min_potential(self) -> None:
        """
        Ensure that the change to potential doesn't cause the potential to exceed the minimum potential.

        :return: None
        """
        self.potential_change = min(self.potential_change, self.minimum_potential)

    def update_potential(self) -> None:
        """
        Update the potential based on any changes calculated.

        :return: None
        """
        self.potential += self.potential_change
        self.potential_change = 0

    def morphic_tansform(self, minimum_potential: float, equilibrium_potential: float, equilibrium_gradient: float,
                         bioelectric_events: Iterable[BioelectricEvent], ligand_sites: Iterable[LigandSite],
                         status: str):
        """
        This method if for applying a morphological change to the cell.

        :param minimum_potential: the new minimum potential
        :param equilibrium_potential: the new equilibrium potential
        :param equilibrium_gradient: The new equilibrium gradient
        :param bioelectric_events: The new bioelectric events
        :param ligand_sites: The new  ligand sites available on the cell
        :param status: The new status of the cell
        :return: None
        """
        self.type = status
        self.equilibrium_potential = equilibrium_potential
        self.equilibrium_gradient = equilibrium_gradient
        self.minimum_potential = minimum_potential
        self.bioelectric_events = bioelectric_events
        self.ligand_sites = {site.ligand: site for site in ligand_sites}
