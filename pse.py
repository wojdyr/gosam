# this file is part of gosam (generator of simple atomistic models)
# Licence (of this file): public domain

"""\
Periodic System of Elements (Z - symbol - name - mass)

Based on:
CRC Handbook of Chemistry & Physics, 63rd edition, 1982-1983
CRC Handbook of Chemistry & Physics, 70th edition, 1989-1990

You may find the same data in cctbx (C++) and atominfo (C).
"""



class Pse:
    def __init__(self, Z, Symbol, Name, Mass):
        self.Z = Z
        self.Symbol = Symbol
        self.Name = Name
        self.Mass = Mass
    def __str__(self):
        return "Z=%s %s (%s)  mass=%s" % (self.Z, self.Name,
                                          self.Symbol, self.Mass)


pse = [
    Pse(   1, "H",  "hydrogen",       1.008 ),
    Pse(   1, "D",  "deuterium",      2.000 ),
    Pse(   2, "He", "helium",         4.003 ),
    Pse(   3, "Li", "lithium",        6.941 ),
    Pse(   4, "Be", "beryllium",      9.012 ),
    Pse(   5, "B",  "boron",         10.811 ),
    Pse(   6, "C",  "carbon",        12.011 ),
    Pse(   7, "N",  "nitrogen",      14.007 ),
    Pse(   8, "O",  "oxygen",        15.999 ),
    Pse(   9, "F",  "fluorine",      18.998 ),
    Pse(  10, "Ne", "neon",          20.180 ),
    Pse(  11, "Na", "sodium",        22.990 ),
    Pse(  12, "Mg", "magnesium",     24.305 ),
    Pse(  13, "Al", "aluminium",     26.982 ),
    Pse(  14, "Si", "silicon",       28.086 ),
    Pse(  15, "P",  "phosphorus",    30.974 ),
    Pse(  16, "S",  "sulphur",       32.066 ),
    Pse(  17, "Cl", "chlorine",      35.452 ),
    Pse(  18, "Ar", "argon",         39.948 ),
    Pse(  19, "K",  "potassium",     39.098 ),
    Pse(  20, "Ca", "calcium",       40.078 ),
    Pse(  21, "Sc", "scandium",      44.956 ),
    Pse(  22, "Ti", "titanium",      47.883 ),
    Pse(  23, "V",  "vanadium",      50.941 ),
    Pse(  24, "Cr", "chromium",      51.996 ),
    Pse(  25, "Mn", "manganese",     54.938 ),
    Pse(  26, "Fe", "iron",          55.847 ),
    Pse(  27, "Co", "cobalt",        58.933 ),
    Pse(  28, "Ni", "nickel",        58.691 ),
    Pse(  29, "Cu", "copper",        63.546 ),
    Pse(  30, "Zn", "zinc",          65.392 ),
    Pse(  31, "Ga", "gallium",       69.723 ),
    Pse(  32, "Ge", "germanium",     72.612 ),
    Pse(  33, "As", "arsenic",       74.922 ),
    Pse(  34, "Se", "selenium",      78.963 ),
    Pse(  35, "Br", "bromine",       79.904 ),
    Pse(  36, "Kr", "krypton",       83.801 ),
    Pse(  37, "Rb", "rubidium",      85.468 ),
    Pse(  38, "Sr", "strontium",     87.621 ),
    Pse(  39, "Y",  "yttrium",       88.906 ),
    Pse(  40, "Zr", "zirconium",     91.224 ),
    Pse(  41, "Nb", "niobium",       92.906 ),
    Pse(  42, "Mo", "molybdenum",    95.941 ),
    Pse(  43, "Tc", "technetium",    98.000 ),
    Pse(  44, "Ru", "ruthenium",    101.072 ),
    Pse(  45, "Rh", "rhodium",      102.905 ),
    Pse(  46, "Pd", "palladium",    106.421 ),
    Pse(  47, "Ag", "silver",       107.868 ),
    Pse(  48, "Cd", "cadmium",      112.411 ),
    Pse(  49, "In", "indium",       114.821 ),
    Pse(  50, "Sn", "tin",          118.710 ),
    Pse(  51, "Sb", "antimony",     121.753 ),
    Pse(  52, "Te", "tellurium",    127.603 ),
    Pse(  53, "I",  "iodine",       126.904 ),
    Pse(  54, "Xe", "xenon",        131.292 ),
    Pse(  55, "Cs", "caesium",      132.905 ),
    Pse(  56, "Ba", "barium",       137.327 ),
    Pse(  57, "La", "lanthanum",    138.906 ),
    Pse(  58, "Ce", "cerium",       140.115 ),
    Pse(  59, "Pr", "praseodymium", 140.908 ),
    Pse(  60, "Nd", "neodymium",    144.243 ),
    Pse(  61, "Pm", "promethium",   145.000 ),
    Pse(  62, "Sm", "samarium",     150.363 ),
    Pse(  63, "Eu", "europium",     151.965 ),
    Pse(  64, "Gd", "gadolinium",   157.253 ),
    Pse(  65, "Tb", "terbium",      158.925 ),
    Pse(  66, "Dy", "dysprosium",   162.503 ),
    Pse(  67, "Ho", "holmium",      164.930 ),
    Pse(  68, "Er", "erbium",       167.263 ),
    Pse(  69, "Tm", "thulium",      168.934 ),
    Pse(  70, "Yb", "ytterbium",    173.043 ),
    Pse(  71, "Lu", "lutetium",     174.967 ),
    Pse(  72, "Hf", "hafnium",      178.492 ),
    Pse(  73, "Ta", "tantalum",     180.948 ),
    Pse(  74, "W",  "tungsten",     183.853 ),
    Pse(  75, "Re", "rhenium",      186.207 ),
    Pse(  76, "Os", "osmium",       190.210 ),
    Pse(  77, "Ir", "iridium",      192.223 ),
    Pse(  78, "Pt", "platinum",     195.083 ),
    Pse(  79, "Au", "gold",         196.967 ),
    Pse(  80, "Hg", "mercury",      200.593 ),
    Pse(  81, "Tl", "thallium",     204.383 ),
    Pse(  82, "Pb", "lead",         207.210 ),
    Pse(  83, "Bi", "bismuth",      208.980 ),
    Pse(  84, "Po", "polonium",     209.000 ),
    Pse(  85, "At", "astatine",     210.000 ),
    Pse(  86, "Rn", "radon",        222.000 ),
    Pse(  87, "Fr", "francium",     223.000 ),
    Pse(  88, "Ra", "radium",       226.025 ),
    Pse(  89, "Ac", "actinium",     227.028 ),
    Pse(  90, "Th", "thorium",      232.038 ),
    Pse(  91, "Pa", "protactinium", 231.035 ),
    Pse(  92, "U",  "uranium",      238.028 ),
    Pse(  93, "Np", "neptunium",    237.048 ),
    Pse(  94, "Pu", "plutonium",    244.000 ),
    Pse(  95, "Am", "americium",    243.000 ),
    Pse(  96, "Cm", "curium",       247.000 ),
    Pse(  97, "Bk", "berkelium",    247.000 ),
    Pse(  98, "Cf", "californium",  251.000 ),
    Pse(  99, "Es", "einsteinium",  254.000 ),
    Pse( 100, "Fm", "fermium",      257.000 ),
    Pse( 101, "Md", "mendelevium",  258.000 ),
    Pse( 102, "No", "nobelium",     259.000 ),
    Pse( 103, "Lr", "lawrencium",   260.000 ),
  ]


# dictionary for easier access using symbols
pse_dict = {}
for i in pse:
    pse_dict[i.Symbol] = i

def get_atom_mass(symbol):
    if symbol in pse_dict:
        return pse_dict[symbol].Mass
    else:
        return 0


