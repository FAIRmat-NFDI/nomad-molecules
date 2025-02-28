import numpy as np
from ase import Atoms
from ase.io import write

from nomad.normalizing import Normalizer
from nomad.datamodel import EntryArchive
from nomad.datamodel.results import Material, System
from molid import query_pubchem_database


class MoleculeNormalizer(Normalizer):
    """Normalizer for molecular data extraction."""
    def normalize(self, archive: EntryArchive, logger=None) -> list:
        self.logger.info("Starting molecular normalization process.")

        if not archive.run:
            print('archive.run: ', archive.run)
            return
        if not EntryArchive().results:
            print('EntryArchive().results: ', EntryArchive().results)
            return

        try:
            # Get the atoms section from the archive
            atoms_data = archive.run[0].system[0].atoms
            atoms = self.create_atoms_object(atoms_data)
            inchikey, molecule_data = self.query_molecule_database(atoms)
            print('self.generate_topology(archive, inchikey, molecule_data): ', self.generate_topology(archive, inchikey, molecule_data))
            return self.generate_topology(archive, inchikey, molecule_data)
        except Exception as e:
            self.logger.error(f"Error in normalization: {e}", exc_info=True)
            print('Error (molecule):', e)
            return

    def create_atoms_object(self, atoms_data) -> Atoms:
        """Creates an ASE Atoms object from NOMAD atomic data."""
        # Get the positions as a numpy array (in meters) and convert to angstrom
        atomic_positions = np.array(atoms_data.positions)
        atomic_positions_angstrom = atomic_positions * 1e10

        # Get the lattice vectors and convert them from meters to angstrom as well
        lattice_vectors = getattr(atoms_data, "lattice_vectors", None)
        if lattice_vectors is not None:
            lattice_vectors = np.array(lattice_vectors, dtype=float) * 1e10
        else:
            lattice_vectors = None

        print('atoms_data.periodic:', atoms_data.periodic)
        atoms = Atoms(
            symbols=atoms_data.labels,
            positions=atomic_positions_angstrom.astype(float),
            cell=lattice_vectors,
            pbc=atoms_data.periodic
        )
        atoms = self.unwrap_atoms(atoms)
        write('atoms.extxyz', atoms)
        return atoms

    def query_molecule_database(self, atoms: Atoms):
        """Queries the local PubChem database using molid."""
        database_file = '../molecule_identification/OpenBabel/pubchem_data_FULL.db'
        try:
            inchikey, molecule_data = query_pubchem_database(atoms, database_file)
            self.logger.info(f"InChIKey: {inchikey}")
            print(f"InChIKey: {inchikey}")
            return inchikey, molecule_data
        except Exception as e:
            self.logger.error(f"Error querying the local PubChem DB: {e}")
            return None, None

    def generate_topology(self, archive, inchikey, molecule_data) -> list:
        """Returns molecular topology data in a format compatible with NOMAD."""
        if not molecule_data:
            self.logger.warning(f"No molecule data found for InChIKey {inchikey}.")
            return []

        # Ensure the archive has the results/material section available.
        if not archive.results:
            archive.results = EntryArchive().results
        if not archive.results.material:
            archive.results.material = archive.results.m_create(Material)
        # import pdb; pdb.set_trace()
        topology_container = archive.m_xpath('results.material.topology')
        if not topology_container:
            topology_container = archive.results.material.topology

        # Create a new System to represent the molecule topology.
        topology_entry = System(
            method='parser',
            label='molecule',
            building_block='molecule'
        )
        # Optionally, you can set additional properties (e.g. symmetry, material_id, etc.)
        # For example:
        # topology_entry.material_id = <your generated id>

        # Append the new system to the topology container.
        topology_container.append(topology_entry)

        return [topology_entry]

    def unwrap_atoms(self, atoms):
        positions = atoms.get_positions(copy=True)
        cell = atoms.get_cell()
        # Check if the cell is non-singular (i.e. has non-zero volume)
        if abs(np.linalg.det(cell)) < 1e-8:
            # If the cell is singular, skip unwrapping.
            return atoms

        # Use the first atom as the reference.
        ref = positions[0]
        for i in range(1, len(positions)):
            disp = positions[i] - ref
            try:
                # Compute the fractional displacement (minimum image)
                frac_disp = np.linalg.solve(cell, disp)
            except np.linalg.LinAlgError:
                # In the unlikely event the cell becomes singular here, skip adjustment for this atom.
                frac_disp = np.zeros_like(disp)
            # Adjust displacement by subtracting the nearest cell translation.
            disp -= cell.dot(np.round(frac_disp))
            positions[i] = ref + disp

        atoms.set_positions(positions)
        atoms.center()  # Optionally recenter the molecule.
        return atoms


