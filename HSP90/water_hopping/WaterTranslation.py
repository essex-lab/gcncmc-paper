import parmed
from simtk import unit
import mdtraj
import numpy as np
import sys, traceback
import math
import copy
import random
import os

class Move(object):

    """Move provides methods for calculating properties on the
    object 'move' (i.e ligand) being perturbed in the NCMC simulation.
    This is the base Move class.
    Ex.
        from blues.ncmc import Model
        ligand = Model(structure, 'LIG')
        ligand.calculateProperties()
    Attributes
    ----------
    """

    def __init__(self):
        """Initialize the Move object
        Currently empy.
        """

    def initializeSystem(self, system, integrator):
        """If the system or integrator needs to be modified to perform the move
        (ex. adding a force) this method is called during the start
        of the simulation to change the system or integrator to accomodate that.
        Parameters
        ----------
        system : simtk.openmm.System object
            System to be modified.
        integrator : simtk.openmm.Integrator object
            Integrator to be modified.
        Returns
        -------
        system : simtk.openmm.System object
            The modified System object.
        integrator : simtk.openmm.Integrator object
            The modified Integrator object.
        """
        new_sys = system
        new_int = integrator
        return new_sys, new_int

    def beforeMove(self, context):
        """This method is called at the start of the NCMC portion if the
        context needs to be checked or modified before performing the move
        at the halfway point.
        Parameters
        ----------
        context: simtk.openmm.Context object
            Context containing the positions to be moved.
        Returns
        -------
        context: simtk.openmm.Context object
            The same input context, but whose context were changed by this function.
        """
        return context

    def afterMove(self, context):
        """This method is called at the end of the NCMC portion if the
        context needs to be checked or modified before performing the move
        at the halfway point.
        Parameters
        ----------
        context: simtk.openmm.Context object
            Context containing the positions to be moved.
        Returns
        -------
        context: simtk.openmm.Context object
            The same input context, but whose context were changed by this function.
        """
        return context

    def _error(self, context):
        """This method is called if running during NCMC portion results
        in an error. This allows portions of the context, such as the
        context parameters that would not be fixed by just reverting the
        positions/velocities of the context.
        Parameters
        ----------
        context: simtk.openmm.Context object
            Context containing the positions to be moved.
        Returns
        -------
        context: simtk.openmm.Context object
            The same input context, but whose context were changed by this function.
        """

        return context

    def move(self, context):
        return context


class WaterTranslationMove(Move):
    """
    Move that translates a random water within a specified radius of the protein's
    center of mass to another point within that radius
    Parameters
    ----------
    structure:
        topology: parmed.Topology
            ParmEd topology object containing atoms of the system.
        water_name: str, optional, default='WAT'
            Residue name of the waters in the system.
        radius: float*unit compatible with simtk.unit.nanometers, optional, default=2.0*unit.nanometers
            Defines the radius within the protein center of mass to choose a water
            and the radius in which to randomly translate that water.
    """

    def __init__(self, structure, water_name='HOH', radius=2.0*unit.nanometers):
        #initialize self attributes
        self.radius = radius #
        self.water_name = water_name
        self.water_residues = []
        protein_atoms = []
        self.before_ncmc_check = True
        self.traj = mdtraj.load('test.pdb')

        residues = list(structure.topology.residues())
        for res in residues:
            if res.name == water_name:
                water_mol = []
                for atom in res.atoms():
                    water_mol.append(atom.index)
                self.water_residues.append(water_mol)
            # Find that atoms of the object that will be used to help specify
            # an atom that will be used to generate the sampling region
            if res.name == 'F09':
                lig = []
                for atom in res.atoms():
                    lig.append(atom.index)
                protein_atoms.append(lig)
        # Choose the first water in the system to act as a bookkeeping device
        self.atom_indices = self.water_residues[0]
        # Indicate the specific atom that the radius will span from
        self.protein_atoms = protein_atoms[0][9]

    def _random_sphere_point(self, radius, origin):
        """
        Function to generate a uniform random point in a sphere of a specified radius.
        Used later on to randomly translate the water molecule.

        Parameters
        ----------
        radius: float
            Defines the radius of the sphere in which a point
            will be uniformly randomly generated.
        """
        # r = specified radius * cubed root of a random number between 0.00 and 0.99999
        r = radius * (np.random.random()**(1./3.))
        # restriction of phi (or azimuth angle) is set from 0 to 2pi.
        # random.uniform allows the values to be chosen w/ an equal probability
        phi = np.random.uniform(0, 2 * np.pi)
        # restriction set from -1 to 1
        costheta = np.random.uniform(-1, 1)
        # calculate theta, the angle between r and Z axis
        theta = np.arccos(costheta)
        # x,y,and z are cartesian coordinates
        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)
        sphere_point = np.array([x, y, z]) * r + origin
        return sphere_point

#    def getMasses(self, topology):
#        """Returns a list of masses of the specified ligand atoms.
#        Parameters
#        ----------
#        topology: parmed.Topology
#            ParmEd topology object containing atoms of the system.
#        """
#        masses = unit.Quantity(np.zeros([int(topology.getNumAtoms()), 1], np.float32), unit.dalton)
#        for idx,atom in enumerate(topology.atoms()):
#            masses[idx] = atom.element._mass #gets the mass of the atom, adds to list (along with index)
#        return masses

    def beforeMove(self, nc_context):
        """
        Temporary function (until multiple alchemical regions are supported),
        which is performed at the beginning of a ncmc iteration. Selects
        a random water within self.radius of the protein's center of mass
        and switches the positions and velocities with the alchemical water
        defined by self.atom_indices, effecitvely duplicating mulitple
        alchemical region support.
        Parameters
        ----------
        nc_context: simtk.openmm Context object
            The context which corresponds to the NCMC simulation.
        """
        start_state = nc_context.getState(getPositions=True, getVelocities=True)
        start_pos = start_state.getPositions(asNumpy=True) #gets starting positions
        start_vel = start_state.getVelocities(asNumpy=True) #gets starting velocities
        switch_pos = np.copy(start_pos)*start_pos.unit
        switch_vel = np.copy(start_vel)*start_vel.unit
        is_inside_sphere = False
        while not is_inside_sphere:
            water_index = np.random.randint(1,len(self.water_residues)) #skip over first water in the system 
            water_choice = self.water_residues[water_index]
            self.traj.xyz[0,:,:] = start_pos;
            pairs = self.traj.topology.select_pairs(np.array(water_choice[0]).flatten(), np.array(self.protein_atoms).flatten())
            water_distance = mdtraj.compute_distances(self.traj, pairs, periodic=True)
            water_dist = np.linalg.norm(water_distance)
            if water_dist <= (self.radius.value_in_unit(unit.nanometers)):
                is_inside_sphere = True
            else:
                pass
        #replace chosen water's positions/velocities with alchemical water
        #the alchemical water is used for bookkeeping
        for i in range(3):
            #set indices of the alchemical waters atoms equal to the indices of the starting positions of the random waters atoms
            switch_pos[self.atom_indices[i]] = start_pos[water_choice[i]]
            #do the same for velocity
            switch_vel[self.atom_indices[i]] = start_vel[water_choice[i]]
            #set indices of the randomly chosen waters atom equal to alchemical waters atom indices. Same w/ velocity
            switch_pos[water_choice[i]] = start_pos[self.atom_indices[i]]
            switch_vel[water_choice[i]] = start_vel[self.atom_indices[i]]
        nc_context.setPositions(switch_pos)
        nc_context.setVelocities(switch_vel)
        return nc_context

    def move(self, context):
        """
        This function is called by the blues.MoveEngine object during a simulation.
        Translates a chosen water to a new randow location within a sphere of self.radius.
        """
        #get the position of the system from the context
        before_move_pos = context.getState(getPositions=True).getPositions(asNumpy=True)
        protein_pos = before_move_pos[self.protein_atoms] #gets the positions from the indices of the atoms in the protein residues in relation to the system
        origin = before_move_pos[self.protein_atoms]
        #Generate uniform random point in a sphere of a specified radius
        point_within_sphere = self._random_sphere_point(self.radius, origin)
        movePos = np.copy(before_move_pos)*before_move_pos.unit #makes a copy of the position of the system from the context
        self.traj.xyz[0,:,:] = movePos;
        pairs = self.traj.topology.select_pairs(np.array(self.atom_indices[0]).flatten(), np.array(self.protein_atoms).flatten())
        water_distance = mdtraj.compute_distances(self.traj, pairs, periodic=True)
        # Replace the oxygens coordinates to the new random location
        movePos[self.atom_indices[0]] = point_within_sphere
        # Get vector distance of the two hydrogens from the oxygen
        H1V = movePos[self.atom_indices[1]] - before_move_pos[self.atom_indices[0]]
        H2V = movePos[self.atom_indices[2]] - before_move_pos[self.atom_indices[0]]
        # Move the hydrogens to new location
        movePos[self.atom_indices[1]] = H1V + point_within_sphere
        movePos[self.atom_indices[2]] = H2V + point_within_sphere
        context.setPositions(movePos)
        return context

    def afterMove(self, nca_context):
        """
        This method is called at the end of the NCMC portion if the
        context needs to be checked or modified before performing the move
        at the halfway point. This function ensures the water is still within
        a sphere of self.radius at the end of the move.
        Parameters
        ----------
        context: simtk.openmm.Context object
            Context containing the positions to be moved.
        Returns
        -------
        context: simtk.openmm.Context object
            The same input context, but whose context were changed by this function.
        """
        before_final_move_pos = nca_context.getState(getPositions=True).getPositions(asNumpy=True)
        # Update positions for distance calculation
        movePos_a = np.copy(before_final_move_pos)*before_final_move_pos.unit
        self.traj.xyz[0,:,:] = movePos_a;
        pairs = self.traj.topology.select_pairs(np.array(self.atom_indices[0]).flatten(), np.array(self.protein_atoms).flatten())
        water_distance = mdtraj.compute_distances(self.traj, pairs, periodic=True)
        # If the water is outside of the radius, reject the move
        if water_distance > np.linalg.norm(self.radius._value):
            nca_context._integrator.setGlobalVariableByName("protocol_work", 999999)
        return nca_context
