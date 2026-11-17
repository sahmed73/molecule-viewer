"""Build the example structures in assets/examples/ with RDKit.

Run from the repository root with RDKit installed (pip install rdkit):

    python3 scripts/make_examples.py

Everyday molecules: SMILES and InChIKeys from PubChem (CIDs below). 3D coordinates
come from RDKit ETKDGv3 conformers optimized with MMFF94; the lowest-energy one is
kept. Each structure is checked against its PubChem InChIKey, so stereochemistry
survives the 3D step.

ZDDP (sec-hexyl): built from a user-supplied SMILES. The two dative S->Zn bonds are
written as single bonds, and zinc is typed as tetrahedral Zn2+ so UFF can optimize
it. The Zn-S geometry is approximate.

Buckminsterfullerene: exact truncated-icosahedron coordinates, since conformer
generators do not handle the closed cage well.

The PAO, ZDDP (dibutyl), and ADPA examples were made separately and are not rebuilt.
"""
import itertools
import math
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

OUT = Path(__file__).resolve().parent.parent / "assets" / "examples"
SEED = 61453

# (file stem, title, PubChem CID, SMILES, InChIKey). Ibuprofen is sold as a racemate
# (CID 3672); the (S) enantiomer is drawn, so its key is the (S) one.
EVERYDAY = [
    ("water", "Water", 962, "O", "XLYOFNOQVPJJNP-UHFFFAOYSA-N"),
    ("carbon-dioxide", "Carbon dioxide", 280, "C(=O)=O", "CURLTUGMZLYLDI-UHFFFAOYSA-N"),
    ("methane", "Methane", 297, "C", "VNWKTOKETHGBQD-UHFFFAOYSA-N"),
    ("ethanol", "Ethanol", 702, "CCO", "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"),
    ("benzene", "Benzene", 241, "C1=CC=CC=C1", "UHOVQNZJYSORNB-UHFFFAOYSA-N"),
    ("caffeine", "Caffeine", 2519, "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "RYYVLZVUVIJVGH-UHFFFAOYSA-N"),
    ("aspirin", "Aspirin", 2244, "CC(=O)OC1=CC=CC=C1C(=O)O", "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"),
    ("paracetamol", "Paracetamol (acetaminophen)", 1983, "CC(=O)NC1=CC=C(C=C1)O", "RZVAJINKPMORJF-UHFFFAOYSA-N"),
    ("ibuprofen", "Ibuprofen, (S) enantiomer", 3672, "CC(C)CC1=CC=C(C=C1)[C@H](C)C(=O)O", None),
    ("nicotine", "Nicotine", 89594, "CN1CCC[C@H]1C2=CN=CC=C2", "SNICXCGAKADSCV-JTQLQIEISA-N"),
    ("dopamine", "Dopamine", 681, "C1=CC(=C(C=C1CCN)O)O", "VYFYYTLLBUKUHU-UHFFFAOYSA-N"),
    ("serotonin", "Serotonin", 5202, "C1=CC2=C(C=C1O)C(=CN2)CCN", "QZAYGJVTTNCVMB-UHFFFAOYSA-N"),
    ("adrenaline", "Adrenaline (epinephrine)", 5816, "CNC[C@@H](C1=CC(=C(C=C1)O)O)O", "UCTWMZQNUQWSLP-VIFPVBQESA-N"),
    ("glucose", "Glucose (beta-D-glucopyranose)", 64689, "C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O", "WQZGKKKJIJFFOK-VFUOTHLCSA-N"),
    ("vitamin-c", "Vitamin C (L-ascorbic acid)", 54670067, "C([C@@H]([C@@H]1C(=C(C(=O)O1)O)O)O)O", "CIWBSHSKHKDKBQ-JLAZNSOCSA-N"),
    ("capsaicin", "Capsaicin", 1548943, "CC(C)/C=C/CCCCC(=O)NCC1=CC(=C(C=C1)O)OC", "YKPUWZUDDOIDPM-SOFGYWHQSA-N"),
    ("cholesterol", "Cholesterol", 5997, "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C", "HVYWMOMLDIMFJA-DPAQBDIFSA-N"),
    ("testosterone", "Testosterone", 6013, "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C", "MUMGGOZAMZWBJJ-DYKIIFRCSA-N"),
    ("penicillin-g", "Penicillin G", 5904, "CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C", "JGSARLDLIJGVTE-MBNYWOFBSA-N"),
    ("morphine", "Morphine", 5288826, "CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O", "BQJCRHHNABKAKU-KBQPJGBKSA-N"),
    ("atp", "ATP (adenosine triphosphate)", 5957, "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N", "ZKHQWZAMYRWXGA-KQYNXXCUSA-N"),
]

ZDDP_SMILES = "[Zn]12(SP(=S->1)(O[C@H](C)CCCC)O[C@H](C)CCCC)SP(=S->2)(O[C@H](C)CCCC)O[C@H](C)CCCC"


def lowest_energy(mol, optimize, conformers=20):
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    params.useRandomCoords = True
    ids = list(AllChem.EmbedMultipleConfs(mol, conformers, params))
    assert ids, "No conformer"
    results = optimize(mol)
    converged = [i for i, (status, _) in enumerate(results) if status == 0] or range(len(ids))
    best = min(converged, key=lambda i: results[i][1])
    out = Chem.Mol(mol)
    keep = Chem.Conformer(mol.GetConformer(ids[best]))
    out.RemoveAllConformers()
    out.AddConformer(keep, assignId=True)
    return out


def centered(mol):
    conf = mol.GetConformer()
    points = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
    cx, cy, cz = (sum(getattr(p, axis) for p in points) / len(points) for axis in "xyz")
    for i, p in enumerate(points):
        conf.SetAtomPosition(i, (p.x - cx, p.y - cy, p.z - cz))
    return mol


def write(stem, title, mol, kekulize=True):
    mol.SetProp("_Name", title)
    OUT.mkdir(parents=True, exist_ok=True)
    block = Chem.MolToMolBlock(centered(mol), kekulize=kekulize)
    assert "V2000" in block
    (OUT / f"{stem}.mol").write_text(block)
    print(f"{stem}.mol: {mol.GetNumAtoms()} atoms")


def everyday(stem, title, cid, smiles, inchikey):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    mol = lowest_energy(mol, lambda m: AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=5000))
    check = Chem.Mol(mol)
    Chem.AssignStereochemistryFrom3D(check)
    expected = inchikey or Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))
    assert Chem.MolToInchiKey(check) == expected, f"{stem}: stereochemistry changed"
    write(stem, f"{title}. PubChem CID {cid}; 3D from RDKit ETKDG + MMFF94", mol)


def zddp_sec_hexyl():
    mol = Chem.RWMol(Chem.AddHs(Chem.MolFromSmiles(ZDDP_SMILES)))
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DATIVE:
            bond.SetBondType(Chem.BondType.SINGLE)
    mol = mol.GetMol()
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
    zinc = next(atom for atom in mol.GetAtoms() if atom.GetSymbol() == "Zn")
    zinc.SetHybridization(Chem.HybridizationType.SP3)  # UFF type Zn3+2
    assert AllChem.UFFHasAllMoleculeParams(mol)
    mol = lowest_energy(mol, lambda m: AllChem.UFFOptimizeMoleculeConfs(m, maxIters=5000), conformers=40)
    conf = mol.GetConformer()
    for neighbor in zinc.GetNeighbors():
        distance = (conf.GetAtomPosition(zinc.GetIdx()) - conf.GetAtomPosition(neighbor.GetIdx())).Length()
        assert 2.0 < distance < 2.6, f"Zn-S distance {distance:.2f}"
    write("zddp-sec-hexyl", "ZDDP: zinc bis[O,O-di(hexan-2-yl) dithiophosphate], (R) carbons; RDKit ETKDG + UFF, approximate Zn-S geometry", mol, kekulize=False)


def fullerene():
    phi = (1 + 5 ** 0.5) / 2
    base = [(0, 1, 3 * phi), (1, 2 + phi, 2 * phi), (phi, 2, phi ** 3)]
    points = set()
    for x, y, z in base:
        for p in [(x, y, z), (y, z, x), (z, x, y)]:  # even permutations
            for signs in itertools.product((1, -1), repeat=3):
                points.add(tuple(round(v * s, 9) + 0.0 for v, s in zip(p, signs)))
    points = sorted(points)
    assert len(points) == 60
    scale = 1.42 / 2  # edges are 2 in these coordinates; C60 bonds average 1.42 Angstrom
    mol = Chem.RWMol()
    for _ in points:
        mol.AddAtom(Chem.Atom(6))
    edges = [(i, j) for i, j in itertools.combinations(range(60), 2) if math.dist(points[i], points[j]) < 2.1]
    assert len(edges) == 90
    for i, j in edges:
        mol.AddBond(i, j, Chem.BondType.SINGLE)
    mol.UpdatePropertyCache(strict=False)
    pentagons = [set(ring) for ring in Chem.GetSymmSSSR(mol) if len(ring) == 5]
    assert len(pentagons) == 12
    for bond in mol.GetBonds():  # Kekule structure: the 30 bonds shared by two hexagons are double
        if not any({bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()} <= ring for ring in pentagons):
            bond.SetBondType(Chem.BondType.DOUBLE)
    mol = mol.GetMol()
    assert sum(b.GetBondType() == Chem.BondType.DOUBLE for b in mol.GetBonds()) == 30
    conf = Chem.Conformer(60)
    for i, p in enumerate(points):
        conf.SetAtomPosition(i, tuple(v * scale for v in p))
    conf.Set3D(True)
    mol.AddConformer(conf)
    write("buckminsterfullerene", "Buckminsterfullerene (C60). PubChem CID 123591; ideal truncated icosahedron, 1.42 Angstrom bonds", mol, kekulize=False)


if __name__ == "__main__":
    for example in EVERYDAY:
        everyday(*example)
    zddp_sec_hexyl()
    fullerene()
