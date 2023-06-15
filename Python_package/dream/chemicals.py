from functools import lru_cache
from itertools import combinations
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS

def pairwise_fingerprint_similarity(smiles, radius=4, n_bits=2048):
    """
    Compute pairwise Tanimoto similarity between fingerprints generated from
    SMILES strings.

    Parameters:
    -----------
    smiles: list of str
        A list of SMILES strings to generate fingerprints from.
    radius: int, optional (default=4)
        The ECFP fingerprint radius.
    n_bits: int, optional (default=2048)
        The length of the ECFP fingerprint in bits.

    Returns:
    --------
    list of tuples
        A list of tuples containing three elements: the first SMILES string,
        the second SMILES string, and the Tanimoto similarity between their
        fingerprints.

    Example:
    --------
    >>> from dream.chemicals import pairwise_fingerprint_similarity
    >>> smiles = ['CCO', 'CCN', 'CCC']
    >>> pairwise_fingerprint_similarity(smiles)
    [('CCO', 'CCN', 0.3333333333333333),
     ('CCO', 'CCC', 0.42857142857142855),
     ('CCN', 'CCC', 0.42857142857142855)]
    """
    fps = {s: ECFP(s, radius=radius, n_bits=n_bits) for s in smiles}
    similarities = [
        (
            smiles1,
            smiles2,
            DataStructs.TanimotoSimilarity(fp1, fp2)
        ) 
        for ((smiles1, fp1), (smiles2, fp2)) in combinations(fps.items(), 2)
    ]
    return similarities

def pairwise_mcs_similarity(smiles, only_heavy_atoms=False):
    """
    Computes the maximum common substructure (MCS) similarity between each pair
    of molecules represented by the input SMILES strings.

    Parameters:
    -----------
        smiles: list of str
            A list of SMILES strings representing the molecules.
        only_heavy_atoms: bool, optional (default=False)
            Consider only heavy atoms (i.e., atoms other than hydrogen) when
            computing the MCS

    Returns:
    --------
        list of tuples
            A list of tuples, where each tuple contains two SMILES strings and
            the corresponding MCS similarity score between the two molecules.

    Example:
    --------
        >>> from dream.chemicals import pairwise_mcs_similarity
        >>> smiles = ['CCO', 'CNC', 'CO']
        >>> pairwise_mcs_similarity(smiles)
        [('CCO', 'CNC', 0.2),
         ('CCO', 'CO', 0.6666666666666666), 
         ('CNC', 'CO', 0.25)]
    """
    mols = {s: mol_from_smiles(s) for s in smiles}
    similarities = [
        (
            smiles1,
            smiles2,
            mcs_overlap(mol1, mol2, only_heavy_atoms=True)
        ) 
        for ((smiles1, mol1), (smiles2, mol2)) in combinations(mols.items(), 2)
    ]
    return similarities

def mcs_overlap(query_mol, target_mol, only_heavy_atoms=False):
    """
    Computes the Maximum Common Substructure (MCS) overlap between two input
    molecules, represented by RDKit Mol objects.

    Parameters:
    -----------
        query_mol: rdkit.Chem.rdchem.Mol
            The query molecule.
        target_mol: rdkit.Chem.rdchem.Mol
            The target molecule.
        only_heavy_atoms: bool, optional (default: False) 
            Consider only heavy atoms (i.e., atoms other than hydrogen) when
            computing the MCS

    Returns:
    --------
        float
            The MCS overlap score between the two molecules, normalized by the
            number of atoms in both molecules.

    Example:
    --------
        >>> from rdkit import Chem
        >>> from dream.chemicals import mcs_overlap
        >>> query_mol = Chem.MolFromSmiles('CCO')
        >>> target_mol = Chem.MolFromSmiles('CNC')
        >>> mcs_overlap(query_mol, target_mol)
        0.2
    """
    if only_heavy_atoms:
        overlap = maximum_common_substructure(
            (query_mol, target_mol), as_mol=True
        ).GetNumHeavyAtoms()
        n_query_atoms = query_mol.GetNumHeavyAtoms()
        n_target_atoms = target_mol.GetNumHeavyAtoms()
    else:
        overlap = maximum_common_substructure((query_mol, target_mol)).numAtoms
        n_query_atoms = query_mol.GetNumAtoms()
        n_target_atoms = target_mol.GetNumAtoms()

    return overlap / (n_query_atoms + n_target_atoms - overlap)

@lru_cache(maxsize=2048)
def mol_from_smiles(smiles):
    """
    Convert a SMILES string to a RDKit Mol molecule object.
    
    Parameters:
    -----------
    smiles : str
        A string representing the SMILES notation of a molecule.

    Returns:
    --------
    rdkit.Chem.Mol
        A RDKit molecule object.
    
    Notes:
    ------
    This function uses the LRU cache to store up to 2048 previously calculated results
    to improve performance when called with the same SMILES string multiple times.
    """
    m = Chem.MolFromSmiles(smiles)
    return m

@lru_cache(maxsize=2048)
def ECFP(smiles, radius=4, n_bits=2048):
    """
    Calculates the Extended-Connectivity Fingerprints (ECFP) of a given SMILES
    string.

    Parameters:
    -----------
    smiles: str
        The SMILES string representing the molecule.

    radius: int, optional (default=4)
        The radius parameter for the ECFP calculation.

    n_bits: int, optional (default=2048)
        The length of the bit vector that represents the ECFP.

    Returns:
    --------
    rdkit.DataStructs.cDataStructs.ExplicitBitVect
        The ECFP fingerprint of the input molecule, as an instance of RDKit's
        ExplicitBitVect class.

    Notes:
    ------
    The function is decorated with the `lru_cache` decorator, which caches the
    results of previous calls for up to 2048 different inputs. This can speed
    up repeated calls to the function with the same input.
    """
    info = {}
    mol = mol_from_smiles(smiles)
    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol, radius=radius, nBits=n_bits, bitInfo=info
    )
    return fp

@lru_cache(maxsize=2048)
def maximum_common_substructure(mols, timeout=1, as_mol=False):
    """
    Find the maximum common substructure of a set of molecules.

    This function uses the RDKit implementation of the maximum common
    substructure algorithm to find the largest common substructure present in a
    set of molecules.

    Parameters:
    -----------
    mols : list of rdkit.Chem.rdchem.Mol
        The list of molecules to compare.
    timeout : float, optional (default=1)
        The maximum amount of time (in seconds) to spend searching for the MCS.
    as_mol : bool, optional (default=False)
        Return the MCS as an RDKit molecule object or as an RDKit MCSResult object.

    Returns:
    --------
    Union[rdkit.Chem.rdchem.Mol, rdkit.Chem.MCSResult]
        If `as_mol` is True, the function returns an RDKit molecule object
        representing the MCS. Otherwise, it returns an RDKit MCSResult object
        with information about the MCS.
    
    Notes:
    ------
    This function is decorated with the `lru_cache` decorator to cache the
    results of previous calls for up to 2048 different inputs and improve
    performance. If the function is called multiple times with the same
    arguments, the cached result will be returned instead of re-running the
    function.
    """
    match = rdFMCS.FindMCS(
        mols,
        timeout=timeout,
        matchValences=True,
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
    )
    if as_mol:
        return Chem.MolFromSmarts(match.smartsString)
    else:
        return match


