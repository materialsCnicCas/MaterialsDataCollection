class PTable:
    def __init__(self):
        self.atom_data = [
            [0, "X", "X", None],  # 0
            [1, "H", "Hydrogen", 1.00794],  # 1
            [2, "He", "Helium", 4.002602],  # 2
            [3, "Li", "Lithium", 6.941],  # 3
            [4, "Be", "Beryllium", 9.012182],  # 4
            [5, "B", "Boron", 10.811],  # 5
            [6, "C", "Carbon", 12.0107],  # 6
            [7, "N", "Nitrogen", 14.0067],  # 7
            [8, "O", "Oxygen", 15.9994],  # 8
            [9, "F", "Fluorine", 18.9984032],  # 9
            [10, "Ne", "Neon", 20.1797],  # 10
            [11, "Na", "Sodium", 22.98976928],  # 11
            [12, "Mg", "Magnesium", 24.3050],  # 12
            [13, "Al", "Aluminium", 26.9815386],  # 13
            [14, "Si", "Silicon", 28.0855],  # 14
            [15, "P", "Phosphorus", 30.973762],  # 15
            [16, "S", "Sulfur", 32.065],  # 16
            [17, "Cl", "Chlorine", 35.453],  # 17
            [18, "Ar", "Argon", 39.948],  # 18
            [19, "K", "Potassium", 39.0983],  # 19
            [20, "Ca", "Calcium", 40.078],  # 20
            [21, "Sc", "Scandium", 44.955912],  # 21
            [22, "Ti", "Titanium", 47.867],  # 22
            [23, "V", "Vanadium", 50.9415],  # 23
            [24, "Cr", "Chromium", 51.9961],  # 24
            [25, "Mn", "Manganese", 54.938045],  # 25
            [26, "Fe", "Iron", 55.845],  # 26
            [27, "Co", "Cobalt", 58.933195],  # 27
            [28, "Ni", "Nickel", 58.6934],  # 28
            [29, "Cu", "Copper", 63.546],  # 29
            [30, "Zn", "Zinc", 65.38],  # 30
            [31, "Ga", "Gallium", 69.723],  # 31
            [32, "Ge", "Germanium", 72.64],  # 32
            [33, "As", "Arsenic", 74.92160],  # 33
            [34, "Se", "Selenium", 78.96],  # 34
            [35, "Br", "Bromine", 79.904],  # 35
            [36, "Kr", "Krypton", 83.798],  # 36
            [37, "Rb", "Rubidium", 85.4678],  # 37
            [38, "Sr", "Strontium", 87.62],  # 38
            [39, "Y", "Yttrium", 88.90585],  # 39
            [40, "Zr", "Zirconium", 91.224],  # 40
            [41, "Nb", "Niobium", 92.90638],  # 41
            [42, "Mo", "Molybdenum", 95.96],  # 42
            [43, "Tc", "Technetium", None],  # 43
            [44, "Ru", "Ruthenium", 101.07],  # 44
            [45, "Rh", "Rhodium", 102.90550],  # 45
            [46, "Pd", "Palladium", 106.42],  # 46
            [47, "Ag", "Silver", 107.8682],  # 47
            [48, "Cd", "Cadmium", 112.411],  # 48
            [49, "In", "Indium", 114.818],  # 49
            [50, "Sn", "Tin", 118.710],  # 50
            [51, "Sb", "Antimony", 121.760],  # 51
            [52, "Te", "Tellurium", 127.60],  # 52
            [53, "I", "Iodine", 126.90447],  # 53
            [54, "Xe", "Xenon", 131.293],  # 54
            [55, "Cs", "Caesium", 132.9054519],  # 55
            [56, "Ba", "Barium", 137.327],  # 56
            [57, "La", "Lanthanum", 138.90547],  # 57
            [58, "Ce", "Cerium", 140.116],  # 58
            [59, "Pr", "Praseodymium", 140.90765],  # 59
            [60, "Nd", "Neodymium", 144.242],  # 60
            [61, "Pm", "Promethium", None],  # 61
            [62, "Sm", "Samarium", 150.36],  # 62
            [63, "Eu", "Europium", 151.964],  # 63
            [64, "Gd", "Gadolinium", 157.25],  # 64
            [65, "Tb", "Terbium", 158.92535],  # 65
            [66, "Dy", "Dysprosium", 162.500],  # 66
            [67, "Ho", "Holmium", 164.93032],  # 67
            [68, "Er", "Erbium", 167.259],  # 68
            [69, "Tm", "Thulium", 168.93421],  # 69
            [70, "Yb", "Ytterbium", 173.054],  # 70
            [71, "Lu", "Lutetium", 174.9668],  # 71
            [72, "Hf", "Hafnium", 178.49],  # 72
            [73, "Ta", "Tantalum", 180.94788],  # 73
            [74, "W", "Tungsten", 183.84],  # 74
            [75, "Re", "Rhenium", 186.207],  # 75
            [76, "Os", "Osmium", 190.23],  # 76
            [77, "Ir", "Iridium", 192.217],  # 77
            [78, "Pt", "Platinum", 195.084],  # 78
            [79, "Au", "Gold", 196.966569],  # 79
            [80, "Hg", "Mercury", 200.59],  # 80
            [81, "Tl", "Thallium", 204.3833],  # 81
            [82, "Pb", "Lead", 207.2],  # 82
            [83, "Bi", "Bismuth", 208.98040],  # 83
            [84, "Po", "Polonium", None],  # 84
            [85, "At", "Astatine", None],  # 85
            [86, "Rn", "Radon", None],  # 86
            [87, "Fr", "Francium", None],  # 87
            [88, "Ra", "Radium", None],  # 88
            [89, "Ac", "Actinium", None],  # 89
            [90, "Th", "Thorium", 232.03806],  # 90
            [91, "Pa", "Protactinium", 231.03588],  # 91
            [92, "U", "Uranium", 238.02891],  # 92
            [93, "Np", "Neptunium", None],  # 93
            [94, "Pu", "Plutonium", None],  # 94
            [95, "Am", "Americium", None],  # 95
            [96, "Cm", "Curium", None],  # 96
            [97, "Bk", "Berkelium", None],  # 97
            [98, "Cf", "Californium", None],  # 98
            [99, "Es", "Einsteinium", None],  # 99
            [100, "Fm", "Fermium", None],  # 100
            [101, "Md", "Mendelevium", None],  # 101
            [102, "No", "Nobelium", None],  # 102
            [103, "Lr", "Lawrencium", None],  # 103
            [104, "Rf", "Rutherfordium", None],  # 104
            [105, "Db", "Dubnium", None],  # 105
            [106, "Sg", "Seaborgium", None],  # 106
            [107, "Bh", "Bohrium", None],  # 107
            [108, "Hs", "Hassium", None],  # 108
            [109, "Mt", "Meitnerium", None],  # 109
            [110, "Ds", "Darmstadtium", None],  # 110
            [111, "Rg", "Roentgenium", None],  # 111
            [112, "Cn", "Copernicium", None],  # 112
            [113, "Uut", "Ununtrium", None],  # 113
            [114, "Uuq", "Ununquadium", None],  # 114
            [115, "Uup", "Ununpentium", None],  # 115
            [116, "Uuh", "Ununhexium", None],  # 116
            [117, "Uus", "Ununseptium", None],  # 117
            [118, "Uuo", "Ununoctium", None],  # 118
        ]
        self.detailed = {
            "Ac": {
                "Atomic mass": 227.0,
                "Atomic no": 89,
                "Atomic orbitals": {
                    "1s": -3443.110367,
                    "2p": -572.7627,
                    "2s": -592.622878,
                    "3d": -119.541743,
                    "3p": -137.654394,
                    "3s": -147.320716,
                    "4d": -23.57061,
                    "4f": -12.278225,
                    "4p": -31.761846,
                    "4s": -36.15826,
                    "5d": -3.222752,
                    "5p": -6.06511,
                    "5s": -7.713078,
                    "6d": -0.137786,
                    "6p": -0.744524,
                    "6s": -1.19698,
                    "7s": -0.126551
                },
                "Atomic radius": 1.95,
                "Atomic radius calculated": "no data",
                "Boiling point": "3573 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "10070 kg m<sup>-3</sup>",
                "Electrical resistivity": "no data 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Rn].6d<sup>1</sup>.7s<sup>2</sup>",
                "Ionic radii": {
                    "3": 1.26
                },
                "Liquid range": "2250 K",
                "Melting point": "1323 K",
                "Mendeleev no": 48,
                "Mineral hardness": "no data",
                "Molar volume": "22.55 cm<sup>3</sup>",
                "Name": "Actinium",
                "Oxidation states": [
                    3
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.26,
                                "ionic_radius": 1.12
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "12 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.1,
                "Youngs modulus": "no data GPa",
                "Metallic radius": 1.878,
                "iupac_ordering": 32,
                "IUPAC ordering": 32
            },
            "Ag": {
                "Atomic mass": 107.8682,
                "Atomic no": 47,
                "Atomic orbitals": {
                    "1s": -900.324578,
                    "2p": -120.913351,
                    "2s": -129.859807,
                    "3d": -13.367803,
                    "3p": -20.06763,
                    "3s": -23.678437,
                    "4d": -0.298706,
                    "4p": -2.086602,
                    "4s": -3.22309,
                    "5s": -0.157407
                },
                "Atomic radius": 1.6,
                "Atomic radius calculated": 1.65,
                "Boiling point": "2435 K",
                "Brinell hardness": "24.5 MN m<sup>-2</sup>",
                "Bulk modulus": "100 GPa",
                "Coefficient of linear thermal expansion": "18.9 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    1
                ],
                "Critical temperature": "no data K",
                "Density of solid": "10490 kg m<sup>-3</sup>",
                "Electrical resistivity": "1.63 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Kr].4d<sup>10</sup>.5s<sup>1</sup>",
                "ICSD oxidation states": [
                    1,
                    2,
                    3
                ],
                "Ionic radii": {
                    "1": 1.29,
                    "2": 1.08,
                    "3": 0.89
                },
                "Liquid range": "1200.07 K",
                "Melting point": "1234.93 K",
                "Mendeleev no": 71,
                "Mineral hardness": "2.5",
                "Molar volume": "10.27 cm<sup>3</sup>",
                "Name": "Silver",
                "Oxidation states": [
                    1,
                    2,
                    3
                ],
                "Poissons ratio": "0.37",
                "Reflectivity": "97 %",
                "Refractive index": "no data",
                "Rigidity modulus": "30 GPa",
                "Shannon radii": {
                    "1": {
                        "II": {
                            "": {
                                "crystal_radius": 0.81,
                                "ionic_radius": 0.67
                            }
                        },
                        "IV": {
                            "": {
                                "crystal_radius": 1.14,
                                "ionic_radius": 1.0
                            }
                        },
                        "IVSQ": {
                            "": {
                                "crystal_radius": 1.16,
                                "ionic_radius": 1.02
                            }
                        },
                        "V": {
                            "": {
                                "crystal_radius": 1.23,
                                "ionic_radius": 1.09
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 1.29,
                                "ionic_radius": 1.15
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.36,
                                "ionic_radius": 1.22
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.42,
                                "ionic_radius": 1.28
                            }
                        }
                    },
                    "2": {
                        "IVSQ": {
                            "": {
                                "crystal_radius": 0.93,
                                "ionic_radius": 0.79
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 1.08,
                                "ionic_radius": 0.94
                            }
                        }
                    },
                    "3": {
                        "IVSQ": {
                            "": {
                                "crystal_radius": 0.81,
                                "ionic_radius": 0.67
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.89,
                                "ionic_radius": 0.75
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "430 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.72,
                "Velocity of sound": "2600 m s<sup>-1</sup>",
                "Vickers hardness": "251 MN m<sup>-2</sup>",
                "X": 1.93,
                "Youngs modulus": "83 GPa",
                "Metallic radius": 1.445,
                "iupac_ordering": 72,
                "IUPAC ordering": 72
            },
            "Al": {
                "Atomic mass": 26.9815386,
                "Atomic no": 13,
                "Atomic orbitals": {
                    "1s": -55.156044,
                    "2p": -2.564018,
                    "2s": -3.934827,
                    "3p": -0.102545,
                    "3s": -0.286883
                },
                "Atomic radius": 1.25,
                "Atomic radius calculated": 1.18,
                "Boiling point": "2792 K",
                "Brinell hardness": "245 MN m<sup>-2</sup>",
                "Bulk modulus": "76 GPa",
                "Coefficient of linear thermal expansion": "23.1 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "2700 kg m<sup>-3</sup>",
                "Electrical resistivity": "2.7 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ne].3s<sup>2</sup>.3p<sup>1</sup>",
                "ICSD oxidation states": [
                    3
                ],
                "Ionic radii": {
                    "3": 0.675
                },
                "Liquid range": "1858.53 K",
                "Melting point": "933.47 K",
                "Mendeleev no": 80,
                "Mineral hardness": "2.75",
                "Molar volume": "10.00 cm<sup>3</sup>",
                "Name": "Aluminum",
                "Oxidation states": [
                    1,
                    3
                ],
                "Poissons ratio": "0.35",
                "Reflectivity": "71 %",
                "Refractive index": "no data",
                "Rigidity modulus": "26 GPa",
                "Shannon radii": {
                    "3": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.53,
                                "ionic_radius": 0.39
                            }
                        },
                        "V": {
                            "": {
                                "crystal_radius": 0.62,
                                "ionic_radius": 0.48
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.675,
                                "ionic_radius": 0.535
                            }
                        }
                    }
                },
                "Superconduction temperature": "1.175 K",
                "Thermal conductivity": "235 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.84,
                "Velocity of sound": "5100 m s<sup>-1</sup>",
                "Vickers hardness": "167 MN m<sup>-2</sup>",
                "X": 1.61,
                "Youngs modulus": "70 GPa",
                "NMR Quadrupole Moment": {
                    "Al-27": 146.6
                },
                "Metallic radius": 1.43,
                "iupac_ordering": 80,
                "IUPAC ordering": 80
            },
            "Am": {
                "Atomic mass": 243.0,
                "Atomic no": 95,
                "Atomic orbitals": "no data",
                "Atomic radius": 1.75,
                "Atomic radius calculated": "no data",
                "Boiling point": "2880 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "no data kg m<sup>-3</sup>",
                "Electrical resistivity": "no data 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Rn].5f<sup>7</sup>.7s<sup>2</sup>",
                "Ionic radii": {
                    "2": 1.4,
                    "3": 1.115,
                    "4": 0.99
                },
                "Liquid range": "1431 K",
                "Melting point": "1449 K",
                "Mendeleev no": 42,
                "Mineral hardness": "no data",
                "Molar volume": "17.63 cm<sup>3</sup>",
                "Name": "Americium",
                "Oxidation states": [
                    2,
                    3,
                    4,
                    5,
                    6
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "2": {
                        "VII": {
                            "": {
                                "crystal_radius": 1.35,
                                "ionic_radius": 1.21
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.4,
                                "ionic_radius": 1.26
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.45,
                                "ionic_radius": 1.31
                            }
                        }
                    },
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.115,
                                "ionic_radius": 0.975
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.23,
                                "ionic_radius": 1.09
                            }
                        }
                    },
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.99,
                                "ionic_radius": 0.85
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.09,
                                "ionic_radius": 0.95
                            }
                        }
                    }
                },
                "Superconduction temperature": "0.6 K",
                "Thermal conductivity": "10 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.3,
                "Youngs modulus": "no data GPa",
                "Metallic radius": 1.73,
                "iupac_ordering": 26,
                "IUPAC ordering": 26
            },
            "Ar": {
                "Atomic mass": 39.948,
                "Atomic no": 18,
                "Atomic orbitals": {
                    "1s": -113.800134,
                    "2p": -8.443439,
                    "2s": -10.794172,
                    "3p": -0.38233,
                    "3s": -0.883384
                },
                "Atomic radius": 0.71,
                "Atomic radius calculated": 0.71,
                "Boiling point": "87.3 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Critical temperature": "150.8 K",
                "Density of solid": "no data kg m<sup>-3</sup>",
                "Electrical resistivity": "no data 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ne].3s<sup>2</sup>.3p<sup>6</sup>",
                "Liquid range": "3.5 K",
                "Max oxidation state": 0.0,
                "Melting point": "83.8 K",
                "Mendeleev no": 3,
                "Min oxidation state": 0.0,
                "Mineral hardness": "no data",
                "Molar volume": "22.56 cm<sup>3</sup>",
                "Name": "Argon",
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "1.000281",
                "Rigidity modulus": "no data GPa",
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "0.01772 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.88,
                "Velocity of sound": "319 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "Youngs modulus": "no data GPa",
                "Metallic radius": "no data",
                "iupac_ordering": 3,
                "IUPAC ordering": 3
            },
            "As": {
                "Atomic mass": 74.9216,
                "Atomic no": 33,
                "Atomic orbitals": {
                    "1s": -423.336658,
                    "2p": -47.527869,
                    "2s": -53.093086,
                    "3d": -1.542767,
                    "3p": -4.851725,
                    "3s": -6.730755,
                    "4p": -0.197497,
                    "4s": -0.52367
                },
                "Atomic radius": 1.15,
                "Atomic radius calculated": 1.14,
                "Boiling point": "887 K",
                "Brinell hardness": "1440 MN m<sup>-2</sup>",
                "Bulk modulus": "22 GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    -3,
                    3,
                    5
                ],
                "Critical temperature": "1700 K",
                "Density of solid": "5727 kg m<sup>-3</sup>",
                "Electrical resistivity": "33 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ar].3d<sup>10</sup>.4s<sup>2</sup>.4p<sup>3</sup>",
                "ICSD oxidation states": [
                    2,
                    3,
                    5,
                    -2,
                    -3,
                    -1
                ],
                "Ionic radii": {
                    "3": 0.72,
                    "5": 0.6
                },
                "Liquid range": "203 K",
                "Melting point": "1090 K",
                "Mendeleev no": 89,
                "Mineral hardness": "3.5",
                "Molar volume": "12.95 cm<sup>3</sup>",
                "Name": "Arsenic",
                "Oxidation states": [
                    -3,
                    2,
                    3,
                    5
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "1.001552",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.72,
                                "ionic_radius": 0.58
                            }
                        }
                    },
                    "5": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.475,
                                "ionic_radius": 0.335
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.6,
                                "ionic_radius": 0.46
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "50 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.85,
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 2.18,
                "Youngs modulus": "8 GPa",
                "Metallic radius": "no data",
                "iupac_ordering": 89,
                "IUPAC ordering": 89
            },
            "At": {
                "Atomic mass": 210.0,
                "Atomic no": 85,
                "Atomic orbitals": {
                    "1s": -3127.390276,
                    "2p": -513.044243,
                    "2s": -531.81835,
                    "3d": -103.060375,
                    "3p": -119.995013,
                    "3s": -129.035542,
                    "4d": -18.295162,
                    "4f": -8.063483,
                    "4p": -25.778264,
                    "4s": -29.809515,
                    "5d": -1.643758,
                    "5p": -4.027061,
                    "5s": -5.453383,
                    "6p": -0.255453,
                    "6s": -0.560189
                },
                "Atomic radius": "no data",
                "Atomic radius calculated": "no data",
                "Boiling point": "no data K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    -1,
                    1
                ],
                "Critical temperature": "no data K",
                "Density of solid": "no data kg m<sup>-3</sup>",
                "Electrical resistivity": "no data 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>14</sup>.5d<sup>10</sup>.6s<sup>2</sup>.6p<sup>5</sup>",
                "Ionic radii": {
                    "7": 0.76
                },
                "Liquid range": "no data K",
                "Melting point": "575 K",
                "Mendeleev no": 96,
                "Mineral hardness": "no data",
                "Molar volume": "no data cm<sup>3</sup>",
                "Name": "Astatine",
                "Oxidation states": [
                    -1,
                    1,
                    3,
                    5
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "7": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.76,
                                "ionic_radius": 0.62
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "2 (estimate)W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 2.02,
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 2.2,
                "Youngs modulus": "no data GPa",
                "Metallic radius": "no data",
                "iupac_ordering": 98,
                "IUPAC ordering": 98
            },
            "Au": {
                "Atomic mass": 196.966569,
                "Atomic no": 79,
                "Atomic orbitals": {
                    "1s": -2683.508245,
                    "2p": -430.725701,
                    "2s": -447.888973,
                    "3d": -81.511751,
                    "3p": -96.707,
                    "3s": -104.824516,
                    "4d": -12.131815,
                    "4f": -3.486824,
                    "4p": -18.578652,
                    "4s": -22.078357,
                    "5d": -0.304738,
                    "5p": -2.002495,
                    "5s": -3.113936,
                    "6s": -0.162334
                },
                "Atomic radius": 1.35,
                "Atomic radius calculated": 1.74,
                "Boiling point": "3129 K",
                "Brinell hardness": "2450 MN m<sup>-2</sup>",
                "Bulk modulus": "220 GPa",
                "Coefficient of linear thermal expansion": "14.2 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "19300 kg m<sup>-3</sup>",
                "Electrical resistivity": "2.2 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>14</sup>.5d<sup>10</sup>.6s<sup>1</sup>",
                "Ionic radii": {
                    "1": 1.51,
                    "3": 0.99,
                    "5": 0.71
                },
                "Liquid range": "1791.67 K",
                "Melting point": "1337.33 K",
                "Mendeleev no": 70,
                "Mineral hardness": "2.5",
                "Molar volume": "10.21 cm<sup>3</sup>",
                "Name": "Gold",
                "Oxidation states": [
                    -1,
                    1,
                    2,
                    3,
                    5
                ],
                "Poissons ratio": "0.44",
                "Reflectivity": "95 %",
                "Refractive index": "no data",
                "Rigidity modulus": "27 GPa",
                "Shannon radii": {
                    "1": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.51,
                                "ionic_radius": 1.37
                            }
                        }
                    },
                    "3": {
                        "IVSQ": {
                            "": {
                                "crystal_radius": 0.82,
                                "ionic_radius": 0.68
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.99,
                                "ionic_radius": 0.85
                            }
                        }
                    },
                    "5": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.71,
                                "ionic_radius": 0.57
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "320 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.66,
                "Velocity of sound": "1740 m s<sup>-1</sup>",
                "Vickers hardness": "216 MN m<sup>-2</sup>",
                "X": 2.54,
                "Youngs modulus": "78 GPa",
                "Metallic radius": 1.442,
                "iupac_ordering": 71,
                "IUPAC ordering": 71
            },
            "B": {
                "Atomic mass": 10.811,
                "Atomic no": 5,
                "Atomic orbitals": {
                    "1s": -6.564347,
                    "2p": -0.136603,
                    "2s": -0.344701
                },
                "Atomic radius": 0.85,
                "Atomic radius calculated": 0.87,
                "Boiling point": "4200 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "320 GPa",
                "Coefficient of linear thermal expansion": "6 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "2460 kg m<sup>-3</sup>",
                "Electrical resistivity": "&gt; 10<sup>12</sup>10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[He].2s<sup>2</sup>.2p<sup>1</sup>",
                "ICSD oxidation states": [
                    3,
                    -3
                ],
                "Ionic radii": {
                    "3": 0.41
                },
                "Liquid range": "1851 K",
                "Melting point": "2349 K",
                "Mendeleev no": 86,
                "Mineral hardness": "9.3",
                "Molar volume": "4.39 cm<sup>3</sup>",
                "Name": "Boron",
                "Oxidation states": [
                    1,
                    2,
                    3
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "3": {
                        "III": {
                            "": {
                                "crystal_radius": 0.15,
                                "ionic_radius": 0.01
                            }
                        },
                        "IV": {
                            "": {
                                "crystal_radius": 0.25,
                                "ionic_radius": 0.11
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.41,
                                "ionic_radius": 0.27
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "27 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.92,
                "Velocity of sound": "16200 m s<sup>-1</sup>",
                "Vickers hardness": "49000 MN m<sup>-2</sup>",
                "X": 2.04,
                "Youngs modulus": "no data GPa",
                "NMR Quadrupole Moment": {
                    "B-10": 84.59,
                    "B-11": 40.59
                },
                "Metallic radius": "no data",
                "iupac_ordering": 81,
                "IUPAC ordering": 81
            },
            "Ba": {
                "Atomic mass": 137.327,
                "Atomic no": 56,
                "Atomic orbitals": {
                    "1s": -1305.743258,
                    "2p": -189.598483,
                    "2s": -200.844444,
                    "3d": -28.528933,
                    "3p": -37.536931,
                    "3s": -42.359434,
                    "4d": -3.432441,
                    "4p": -6.497622,
                    "4s": -8.257061,
                    "5p": -0.698605,
                    "5s": -1.157159,
                    "6s": -0.118967
                },
                "Atomic radius": 2.15,
                "Atomic radius calculated": 2.53,
                "Boiling point": "2143 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "9.6 GPa",
                "Coefficient of linear thermal expansion": "20.6 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    2
                ],
                "Critical temperature": "no data K",
                "Density of solid": "3510 kg m<sup>-3</sup>",
                "Electrical resistivity": "34 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].6s<sup>2</sup>",
                "ICSD oxidation states": [
                    2
                ],
                "Ionic radii": {
                    "2": 1.49
                },
                "Liquid range": "1143 K",
                "Melting point": "1000 K",
                "Mendeleev no": 14,
                "Mineral hardness": "1.25",
                "Molar volume": "38.16 cm<sup>3</sup>",
                "Name": "Barium",
                "Oxidation states": [
                    2
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "4.9 GPa",
                "Shannon radii": {
                    "2": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.49,
                                "ionic_radius": 1.35
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.52,
                                "ionic_radius": 1.38
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.56,
                                "ionic_radius": 1.42
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.61,
                                "ionic_radius": 1.47
                            }
                        },
                        "X": {
                            "": {
                                "crystal_radius": 1.66,
                                "ionic_radius": 1.52
                            }
                        },
                        "XI": {
                            "": {
                                "crystal_radius": 1.71,
                                "ionic_radius": 1.57
                            }
                        },
                        "XII": {
                            "": {
                                "crystal_radius": 1.75,
                                "ionic_radius": 1.61
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "18 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 2.68,
                "Velocity of sound": "1620 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 0.89,
                "Youngs modulus": "13 GPa",
                "Metallic radius": 2.236,
                "iupac_ordering": 13,
                "IUPAC ordering": 13
            },
            "Be": {
                "Atomic mass": 9.012182,
                "Atomic no": 4,
                "Atomic orbitals": {
                    "1s": -3.856411,
                    "2s": -0.205744
                },
                "Atomic radius": 1.05,
                "Atomic radius calculated": 1.12,
                "Boiling point": "2742 K",
                "Brinell hardness": "600 MN m<sup>-2</sup>",
                "Bulk modulus": "130 GPa",
                "Coefficient of linear thermal expansion": "11.3 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    2
                ],
                "Critical temperature": "no data K",
                "Density of solid": "1848 kg m<sup>-3</sup>",
                "Electrical resistivity": "3.8 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[He].2s<sup>2</sup>",
                "ICSD oxidation states": [
                    2
                ],
                "Ionic radii": {
                    "2": 0.59
                },
                "Liquid range": "1182 K",
                "Melting point": "1560 K",
                "Mendeleev no": 77,
                "Mineral hardness": "5.5",
                "Molar volume": "4.85 cm<sup>3</sup>",
                "Name": "Beryllium",
                "Oxidation states": [
                    2
                ],
                "Poissons ratio": "0.032",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "132 GPa",
                "Shannon radii": {
                    "2": {
                        "III": {
                            "": {
                                "crystal_radius": 0.3,
                                "ionic_radius": 0.16
                            }
                        },
                        "IV": {
                            "": {
                                "crystal_radius": 0.41,
                                "ionic_radius": 0.27
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.59,
                                "ionic_radius": 0.45
                            }
                        }
                    }
                },
                "Superconduction temperature": "0.026 K",
                "Thermal conductivity": "190 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.53,
                "Velocity of sound": "13000 m s<sup>-1</sup>",
                "Vickers hardness": "1670 MN m<sup>-2</sup>",
                "X": 1.57,
                "Youngs modulus": "287 GPa",
                "NMR Quadrupole Moment": {
                    "Be-9": 52.88
                },
                "Metallic radius": 1.12,
                "iupac_ordering": 17,
                "IUPAC ordering": 17
            },
            "Bi": {
                "Atomic mass": 208.9804,
                "Atomic no": 83,
                "Atomic orbitals": {
                    "1s": -2975.550959,
                    "2p": -484.716359,
                    "2s": -502.950758,
                    "3d": -95.532476,
                    "3p": -111.883393,
                    "3s": -120.613998,
                    "4d": -16.084817,
                    "4f": -6.382744,
                    "4p": -23.218641,
                    "4s": -27.07034,
                    "5d": -1.139408,
                    "5p": -3.293637,
                    "5s": -4.611934,
                    "6p": -0.180198,
                    "6s": -0.426129
                },
                "Atomic radius": 1.6,
                "Atomic radius calculated": 1.43,
                "Boiling point": "1837 K",
                "Brinell hardness": "94.2 MN m<sup>-2</sup>",
                "Bulk modulus": "31 GPa",
                "Coefficient of linear thermal expansion": "13.4 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "9780 kg m<sup>-3</sup>",
                "Electrical resistivity": "130 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>14</sup>.5d<sup>10</sup>.6s<sup>2</sup>.6p<sup>3</sup>",
                "ICSD oxidation states": [
                    1,
                    2,
                    3,
                    5
                ],
                "Ionic radii": {
                    "3": 1.17,
                    "5": 0.9
                },
                "Liquid range": "1292.6 K",
                "Melting point": "544.4 K",
                "Mendeleev no": 87,
                "Mineral hardness": "2.25",
                "Molar volume": "21.31 cm<sup>3</sup>",
                "Name": "Bismuth",
                "Oxidation states": [
                    -3,
                    3,
                    5
                ],
                "Poissons ratio": "0.33",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "12 GPa",
                "Shannon radii": {
                    "3": {
                        "V": {
                            "": {
                                "crystal_radius": 1.1,
                                "ionic_radius": 0.96
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 1.17,
                                "ionic_radius": 1.03
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.31,
                                "ionic_radius": 1.17
                            }
                        }
                    },
                    "5": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.9,
                                "ionic_radius": 0.76
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "8 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 2.07,
                "Velocity of sound": "1790 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 2.02,
                "Youngs modulus": "32 GPa",
                "Metallic radius": 1.82,
                "iupac_ordering": 87,
                "IUPAC ordering": 87
            },
            "Bk": {
                "Atomic mass": 247.0,
                "Atomic no": 97,
                "Atomic orbitals": "no data",
                "Atomic radius": "no data",
                "Atomic radius calculated": "no data",
                "Boiling point": "no data K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "14780 kg m<sup>-3</sup>",
                "Electrical resistivity": "no data 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Rn].5f<sup>9</sup>.7s<sup>2</sup>",
                "Ionic radii": {
                    "3": 1.1,
                    "4": 0.97
                },
                "Liquid range": "no data K",
                "Melting point": "1259 K",
                "Mendeleev no": 40,
                "Mineral hardness": "no data",
                "Molar volume": "16.84 cm<sup>3</sup>",
                "Name": "Berkelium",
                "Oxidation states": [
                    3,
                    4
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.1,
                                "ionic_radius": 0.96
                            }
                        }
                    },
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.97,
                                "ionic_radius": 0.83
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.07,
                                "ionic_radius": 0.93
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "10 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.3,
                "Youngs modulus": "no data GPa",
                "Metallic radius": 1.703,
                "iupac_ordering": 24,
                "IUPAC ordering": 24
            },
            "Br": {
                "Atomic mass": 79.904,
                "Atomic no": 35,
                "Atomic orbitals": {
                    "1s": -480.182643,
                    "2p": -55.67796,
                    "2s": -61.710022,
                    "3d": -2.52211,
                    "3p": -6.298805,
                    "3s": -8.409057,
                    "4p": -0.295334,
                    "4s": -0.720066
                },
                "Atomic radius": 1.15,
                "Atomic radius calculated": 0.94,
                "Boiling point": "332 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "1.9 GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    -1,
                    1,
                    3,
                    5,
                    7
                ],
                "Critical temperature": "586 K",
                "Density of solid": "no data kg m<sup>-3</sup>",
                "Electrical resistivity": "&gt; 10<sup>18</sup>10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ar].3d<sup>10</sup>.4s<sup>2</sup>.4p<sup>5</sup>",
                "ICSD oxidation states": [
                    5,
                    -1
                ],
                "Ionic radii": {
                    "-1": 1.82,
                    "3": 0.73,
                    "5": 0.45,
                    "7": 0.53
                },
                "Liquid range": "66.2 K",
                "Melting point": "265.8 K",
                "Mendeleev no": 98,
                "Mineral hardness": "no data",
                "Molar volume": "19.78 cm<sup>3</sup>",
                "Name": "Bromine",
                "Oxidation states": [
                    -1,
                    1,
                    3,
                    4,
                    5,
                    7
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "1.001132",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "-1": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.82,
                                "ionic_radius": 1.96
                            }
                        }
                    },
                    "3": {
                        "IVSQ": {
                            "": {
                                "crystal_radius": 0.73,
                                "ionic_radius": 0.59
                            }
                        }
                    },
                    "5": {
                        "IIIPY": {
                            "": {
                                "crystal_radius": 0.45,
                                "ionic_radius": 0.31
                            }
                        }
                    },
                    "7": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.39,
                                "ionic_radius": 0.25
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.53,
                                "ionic_radius": 0.39
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "0.12 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.85,
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 2.96,
                "Youngs modulus": "no data GPa",
                "Metallic radius": 1.14,
                "iupac_ordering": 100,
                "IUPAC ordering": 100
            },
            "C": {
                "Atomic mass": 12.0107,
                "Atomic no": 6,
                "Atomic orbitals": {
                    "1s": -9.947718,
                    "2p": -0.199186,
                    "2s": -0.500866
                },
                "Atomic radius": 0.7,
                "Atomic radius calculated": 0.67,
                "Boiling point": "4300 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "33 GPa",
                "Coefficient of linear thermal expansion": "7.1 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    -4,
                    4
                ],
                "Critical temperature": "no data K",
                "Density of solid": "2267 kg m<sup>-3</sup>",
                "Electrical resistivity": "about 1000 - direction dependent10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[He].2s<sup>2</sup>.2p<sup>2</sup>",
                "ICSD oxidation states": [
                    2,
                    3,
                    4,
                    -4,
                    -3,
                    -2
                ],
                "Ionic radii": {
                    "4": 0.3
                },
                "Liquid range": "500 K",
                "Melting point": "3800 K",
                "Mendeleev no": 95,
                "Mineral hardness": "0.5 (graphite; diamond is 10.0)(no units)",
                "Molar volume": "5.29 cm<sup>3</sup>",
                "Name": "Carbon",
                "Oxidation states": [
                    -4,
                    -3,
                    -2,
                    -1,
                    1,
                    2,
                    3,
                    4
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "27 %",
                "Refractive index": "2.417 (diamond)(no units)",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "4": {
                        "III": {
                            "": {
                                "crystal_radius": 0.06,
                                "ionic_radius": -0.08
                            }
                        },
                        "IV": {
                            "": {
                                "crystal_radius": 0.29,
                                "ionic_radius": 0.15
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.3,
                                "ionic_radius": 0.16
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "140 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.7,
                "Velocity of sound": "18350 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 2.55,
                "Youngs modulus": "no data GPa",
                "NMR Quadrupole Moment": {
                    "C-11": 33.27
                },
                "Metallic radius": "no data",
                "iupac_ordering": 86,
                "IUPAC ordering": 86
            },
            "Ca": {
                "Atomic mass": 40.078,
                "Atomic no": 20,
                "Atomic orbitals": {
                    "1s": -143.935181,
                    "2p": -12.285376,
                    "2s": -15.046905,
                    "3p": -1.030572,
                    "3s": -1.706331,
                    "4s": -0.141411
                },
                "Atomic radius": 1.8,
                "Atomic radius calculated": 1.94,
                "Boiling point": "1757 K",
                "Brinell hardness": "167 MN m<sup>-2</sup>",
                "Bulk modulus": "17 GPa",
                "Coefficient of linear thermal expansion": "22.3 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    2
                ],
                "Critical temperature": "no data K",
                "Density of solid": "1550 kg m<sup>-3</sup>",
                "Electrical resistivity": "3.4 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ar].4s<sup>2</sup>",
                "ICSD oxidation states": [
                    2
                ],
                "Ionic radii": {
                    "2": 1.14
                },
                "Liquid range": "642 K",
                "Melting point": "1115 K",
                "Mendeleev no": 16,
                "Mineral hardness": "1.75",
                "Molar volume": "26.20 cm<sup>3</sup>",
                "Name": "Calcium",
                "Oxidation states": [
                    2
                ],
                "Poissons ratio": "0.31",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "7.4 GPa",
                "Shannon radii": {
                    "2": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.14,
                                "ionic_radius": 1.0
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.2,
                                "ionic_radius": 1.06
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.26,
                                "ionic_radius": 1.12
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.32,
                                "ionic_radius": 1.18
                            }
                        },
                        "X": {
                            "": {
                                "crystal_radius": 1.37,
                                "ionic_radius": 1.23
                            }
                        },
                        "XII": {
                            "": {
                                "crystal_radius": 1.48,
                                "ionic_radius": 1.34
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "200 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 2.31,
                "Velocity of sound": "3810 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.0,
                "Youngs modulus": "20 GPa",
                "NMR Quadrupole Moment": {
                    "Ca-41": -66.5,
                    "Ca-43": -40.8
                },
                "Metallic radius": 1.976,
                "iupac_ordering": 15,
                "IUPAC ordering": 15
            },
            "Cd": {
                "Atomic mass": 112.411,
                "Atomic no": 48,
                "Atomic orbitals": {
                    "1s": -941.476646,
                    "2p": -127.63512,
                    "2s": -136.83249,
                    "3d": -14.685252,
                    "3p": -21.637522,
                    "3s": -25.379908,
                    "4d": -0.47053,
                    "4p": -2.39526,
                    "4s": -3.596069,
                    "5s": -0.204228
                },
                "Atomic radius": 1.55,
                "Atomic radius calculated": 1.61,
                "Boiling point": "1040 K",
                "Brinell hardness": "203 MN m<sup>-2</sup>",
                "Bulk modulus": "42 GPa",
                "Coefficient of linear thermal expansion": "30.8 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    2
                ],
                "Critical temperature": "no data K",
                "Density of solid": "8650 kg m<sup>-3</sup>",
                "Electrical resistivity": "7 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Kr].4d<sup>10</sup>.5s<sup>2</sup>",
                "ICSD oxidation states": [
                    2
                ],
                "Ionic radii": {
                    "2": 1.09
                },
                "Liquid range": "445.78 K",
                "Melting point": "594.22 K",
                "Mendeleev no": 75,
                "Mineral hardness": "2.0",
                "Molar volume": "13.00 cm<sup>3</sup>",
                "Name": "Cadmium",
                "Oxidation states": [
                    1,
                    2
                ],
                "Poissons ratio": "0.30",
                "Reflectivity": "67 %",
                "Refractive index": "no data",
                "Rigidity modulus": "19 GPa",
                "Shannon radii": {
                    "2": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.92,
                                "ionic_radius": 0.78
                            }
                        },
                        "V": {
                            "": {
                                "crystal_radius": 1.01,
                                "ionic_radius": 0.87
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 1.09,
                                "ionic_radius": 0.95
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.17,
                                "ionic_radius": 1.03
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.24,
                                "ionic_radius": 1.1
                            }
                        },
                        "XII": {
                            "": {
                                "crystal_radius": 1.45,
                                "ionic_radius": 1.31
                            }
                        }
                    }
                },
                "Superconduction temperature": "0.517 K",
                "Thermal conductivity": "97 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.58,
                "Velocity of sound": "2310 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.69,
                "Youngs modulus": "50 GPa",
                "Metallic radius": 1.51,
                "iupac_ordering": 75,
                "IUPAC ordering": 75
            },
            "Ce": {
                "Atomic mass": 140.116,
                "Atomic no": 58,
                "Atomic orbitals": {
                    "1s": -1406.148284,
                    "2p": -206.925148,
                    "2s": -218.684842,
                    "3d": -32.412569,
                    "3p": -41.938282,
                    "3s": -47.035283,
                    "4d": -4.192548,
                    "4f": -0.337442,
                    "4p": -7.532106,
                    "4s": -9.432744,
                    "5d": -0.14055,
                    "5p": -0.85011,
                    "5s": -1.369728,
                    "6s": -0.133974
                },
                "Atomic radius": 1.85,
                "Atomic radius calculated": "no data",
                "Boiling point": "3633 K",
                "Brinell hardness": "412 MN m<sup>-2</sup>",
                "Bulk modulus": "22 GPa",
                "Coefficient of linear thermal expansion": "6.3 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3,
                    4
                ],
                "Critical temperature": "no data K",
                "Density of solid": "6689 kg m<sup>-3</sup>",
                "Electrical resistivity": "74 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>1</sup>.5d<sup>1</sup>.6s<sup>2</sup>",
                "ICSD oxidation states": [
                    3,
                    4
                ],
                "Ionic radii": {
                    "3": 1.15,
                    "4": 1.01
                },
                "Liquid range": "2565 K",
                "Melting point": "1068 K",
                "Mendeleev no": 32,
                "Mineral hardness": "2.5",
                "Molar volume": "20.69 cm<sup>3</sup>",
                "Name": "Cerium",
                "Oxidation states": [
                    2,
                    3,
                    4
                ],
                "Poissons ratio": "0.24",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "14 GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.15,
                                "ionic_radius": 1.01
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.21,
                                "ionic_radius": 1.07
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.283,
                                "ionic_radius": 1.143
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.336,
                                "ionic_radius": 1.196
                            }
                        },
                        "X": {
                            "": {
                                "crystal_radius": 1.39,
                                "ionic_radius": 1.25
                            }
                        },
                        "XII": {
                            "": {
                                "crystal_radius": 1.48,
                                "ionic_radius": 1.34
                            }
                        }
                    },
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.01,
                                "ionic_radius": 0.87
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.11,
                                "ionic_radius": 0.97
                            }
                        },
                        "X": {
                            "": {
                                "crystal_radius": 1.21,
                                "ionic_radius": 1.07
                            }
                        },
                        "XII": {
                            "": {
                                "crystal_radius": 1.28,
                                "ionic_radius": 1.14
                            }
                        }
                    }
                },
                "Superconduction temperature": "0.022  (under pressure)K",
                "Thermal conductivity": "11 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "2100 m s<sup>-1</sup>",
                "Vickers hardness": "270 MN m<sup>-2</sup>",
                "X": 1.12,
                "Youngs modulus": "34 GPa",
                "Metallic radius": 1.707,
                "iupac_ordering": 46,
                "IUPAC ordering": 46
            },
            "Cf": {
                "Atomic mass": 251.0,
                "Atomic no": 98,
                "Atomic orbitals": "no data",
                "Atomic radius": "no data",
                "Atomic radius calculated": "no data",
                "Boiling point": "no data K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "15100 kg m<sup>-3</sup>",
                "Electrical resistivity": "no data 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Rn].5f<sup>10</sup>.7s<sup>2</sup>",
                "Ionic radii": {
                    "3": 1.09,
                    "4": 0.961
                },
                "Liquid range": "no data K",
                "Melting point": "1173 K",
                "Mendeleev no": 39,
                "Mineral hardness": "no data",
                "Molar volume": "16.50 cm<sup>3</sup>",
                "Name": "Californium",
                "Oxidation states": [
                    2,
                    3,
                    4
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.09,
                                "ionic_radius": 0.95
                            }
                        }
                    },
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.961,
                                "ionic_radius": 0.821
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.06,
                                "ionic_radius": 0.92
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "no data W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.3,
                "Youngs modulus": "no data GPa",
                "Metallic radius": 1.86,
                "iupac_ordering": 23,
                "IUPAC ordering": 23
            },
            "Cl": {
                "Atomic mass": 35.453,
                "Atomic no": 17,
                "Atomic orbitals": {
                    "1s": -100.369229,
                    "2p": -7.039982,
                    "2s": -9.187993,
                    "3p": -0.32038,
                    "3s": -0.754458
                },
                "Atomic radius": 1.0,
                "Atomic radius calculated": 0.79,
                "Boiling point": "239.11 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "1.1 (liquid)GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    -1,
                    1,
                    3,
                    5,
                    7
                ],
                "Critical temperature": "417 K",
                "Density of solid": "no data kg m<sup>-3</sup>",
                "Electrical resistivity": "&gt; 10<sup>10</sup>10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ne].3s<sup>2</sup>.3p<sup>5</sup>",
                "ICSD oxidation states": [
                    -1
                ],
                "Ionic radii": {
                    "-1": 1.67,
                    "5": 0.26,
                    "7": 0.41
                },
                "Liquid range": "67.51 K",
                "Melting point": "171.6 K",
                "Mendeleev no": 99,
                "Mineral hardness": "no data",
                "Molar volume": "17.39 cm<sup>3</sup>",
                "Name": "Chlorine",
                "Oxidation states": [
                    -1,
                    1,
                    2,
                    3,
                    4,
                    5,
                    6,
                    7
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "1.000773",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "-1": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.67,
                                "ionic_radius": 1.81
                            }
                        }
                    },
                    "5": {
                        "IIIPY": {
                            "": {
                                "crystal_radius": 0.26,
                                "ionic_radius": 0.12
                            }
                        }
                    },
                    "7": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.22,
                                "ionic_radius": 0.08
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.41,
                                "ionic_radius": 0.27
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "0.0089 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.75,
                "Velocity of sound": "206 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 3.16,
                "Youngs modulus": "no data GPa",
                "NMR Quadrupole Moment": {
                    "Cl-35": -81.65,
                    "Cl-37": -64.35
                },
                "Metallic radius": "no data",
                "iupac_ordering": 101,
                "IUPAC ordering": 101
            },
            "Cm": {
                "Atomic mass": 247.0,
                "Atomic no": 96,
                "Atomic orbitals": "no data",
                "Atomic radius": "no data",
                "Atomic radius calculated": "no data",
                "Boiling point": "3383 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "13510 kg m<sup>-3</sup>",
                "Electrical resistivity": "no data 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Rn].5f<sup>7</sup>.6d<sup>1</sup>.7s<sup>2</sup>",
                "Ionic radii": {
                    "3": 1.11,
                    "4": 0.99
                },
                "Liquid range": "1770 K",
                "Melting point": "1613 K",
                "Mendeleev no": 41,
                "Mineral hardness": "no data",
                "Molar volume": "18.05 cm<sup>3</sup>",
                "Name": "Curium",
                "Oxidation states": [
                    3,
                    4
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.11,
                                "ionic_radius": 0.97
                            }
                        }
                    },
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.99,
                                "ionic_radius": 0.85
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.09,
                                "ionic_radius": 0.95
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "8.8 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.3,
                "Youngs modulus": "no data GPa",
                "Metallic radius": 1.743,
                "iupac_ordering": 25,
                "IUPAC ordering": 25
            },
            "Co": {
                "Atomic mass": 58.933195,
                "Atomic no": 27,
                "Atomic orbitals": {
                    "1s": -275.616639,
                    "2p": -28.152095,
                    "2s": -32.379758,
                    "3d": -0.322368,
                    "3p": -2.388285,
                    "3s": -3.651812,
                    "4s": -0.204497
                },
                "Atomic radius": 1.35,
                "Atomic radius calculated": 1.52,
                "Boiling point": "3200 K",
                "Brinell hardness": "700 MN m<sup>-2</sup>",
                "Bulk modulus": "180 GPa",
                "Coefficient of linear thermal expansion": "13.0 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    2,
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "8900 kg m<sup>-3</sup>",
                "Electrical resistivity": "6 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ar].3d<sup>7</sup>.4s<sup>2</sup>",
                "ICSD oxidation states": [
                    1,
                    2,
                    3,
                    4
                ],
                "Ionic radii": {
                    "2": 0.885,
                    "3": 0.75,
                    "4": 0.67
                },
                "Ionic radii hs": {
                    "2": 0.885,
                    "3": 0.75,
                    "4": 0.67
                },
                "Ionic radii ls": {
                    "2": 0.79,
                    "3": 0.685
                },
                "Liquid range": "1432 K",
                "Melting point": "1768 K",
                "Mendeleev no": 64,
                "Mineral hardness": "5.0",
                "Molar volume": "6.67 cm<sup>3</sup>",
                "Name": "Cobalt",
                "Oxidation states": [
                    -1,
                    1,
                    2,
                    3,
                    4,
                    5
                ],
                "Poissons ratio": "0.31",
                "Reflectivity": "67 %",
                "Refractive index": "no data",
                "Rigidity modulus": "75 GPa",
                "Shannon radii": {
                    "2": {
                        "IV": {
                            "High Spin": {
                                "crystal_radius": 0.72,
                                "ionic_radius": 0.58
                            }
                        },
                        "V": {
                            "": {
                                "crystal_radius": 0.81,
                                "ionic_radius": 0.67
                            }
                        },
                        "VI": {
                            "Low Spin": {
                                "crystal_radius": 0.79,
                                "ionic_radius": 0.65
                            },
                            "High Spin": {
                                "crystal_radius": 0.885,
                                "ionic_radius": 0.745
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.04,
                                "ionic_radius": 0.9
                            }
                        }
                    },
                    "3": {
                        "VI": {
                            "High Spin": {
                                "crystal_radius": 0.75,
                                "ionic_radius": 0.61
                            },
                            "Low Spin": {
                                "crystal_radius": 0.685,
                                "ionic_radius": 0.545
                            }
                        }
                    },
                    "4": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.54,
                                "ionic_radius": 0.4
                            }
                        },
                        "VI": {
                            "High Spin": {
                                "crystal_radius": 0.67,
                                "ionic_radius": 0.53
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "100 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "4720 m s<sup>-1</sup>",
                "Vickers hardness": "1043 MN m<sup>-2</sup>",
                "X": 1.88,
                "Youngs modulus": "209 GPa",
                "NMR Quadrupole Moment": {
                    "Co-59": 420.3
                },
                "Metallic radius": 1.25,
                "iupac_ordering": 67,
                "IUPAC ordering": 67
            },
            "Cr": {
                "Atomic mass": 51.9961,
                "Atomic no": 24,
                "Atomic orbitals": {
                    "1s": -213.881191,
                    "2p": -20.526273,
                    "2s": -24.113457,
                    "3d": -0.118123,
                    "3p": -1.65423,
                    "3s": -2.649085,
                    "4s": -0.150445
                },
                "Atomic radius": 1.4,
                "Atomic radius calculated": 1.66,
                "Boiling point": "2944 K",
                "Brinell hardness": "1120 MN m<sup>-2</sup>",
                "Bulk modulus": "160 GPa",
                "Coefficient of linear thermal expansion": "4.9 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3,
                    6
                ],
                "Critical temperature": "no data K",
                "Density of solid": "7140 kg m<sup>-3</sup>",
                "Electrical resistivity": "12.7 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ar].3d<sup>5</sup>.4s<sup>1</sup>",
                "ICSD oxidation states": [
                    2,
                    3,
                    4,
                    5,
                    6
                ],
                "Ionic radii": {
                    "2": 0.94
                },
                "Ionic radii hs": {
                    "2": 0.94
                },
                "Ionic radii ls": {
                    "2": 0.87,
                    "3": 0.755,
                    "4": 0.69,
                    "5": 0.63,
                    "6": 0.58
                },
                "Liquid range": "764 K",
                "Melting point": "2180 K",
                "Mendeleev no": 57,
                "Mineral hardness": "8.5",
                "Molar volume": "7.23 cm<sup>3</sup>",
                "Name": "Chromium",
                "Oxidation states": [
                    -2,
                    -1,
                    1,
                    2,
                    3,
                    4,
                    5,
                    6
                ],
                "Poissons ratio": "0.21",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "115 GPa",
                "Shannon radii": {
                    "2": {
                        "VI": {
                            "Low Spin": {
                                "crystal_radius": 0.87,
                                "ionic_radius": 0.73
                            },
                            "High Spin": {
                                "crystal_radius": 0.94,
                                "ionic_radius": 0.8
                            }
                        }
                    },
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.755,
                                "ionic_radius": 0.615
                            }
                        }
                    },
                    "4": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.55,
                                "ionic_radius": 0.41
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.69,
                                "ionic_radius": 0.55
                            }
                        }
                    },
                    "5": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.485,
                                "ionic_radius": 0.345
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.63,
                                "ionic_radius": 0.49
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 0.71,
                                "ionic_radius": 0.57
                            }
                        }
                    },
                    "6": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.4,
                                "ionic_radius": 0.26
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.58,
                                "ionic_radius": 0.44
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "94 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "5940 m s<sup>-1</sup>",
                "Vickers hardness": "1060 MN m<sup>-2</sup>",
                "X": 1.66,
                "Youngs modulus": "279 GPa",
                "NMR Quadrupole Moment": {
                    "Cr-53": -150.5
                },
                "Metallic radius": 1.285,
                "iupac_ordering": 58,
                "IUPAC ordering": 58
            },
            "Cs": {
                "Atomic mass": 132.9054519,
                "Atomic no": 55,
                "Atomic orbitals": {
                    "1s": -1256.738791,
                    "2p": -180.995344,
                    "2s": -191.981873,
                    "3d": -26.418398,
                    "3p": -35.166423,
                    "3s": -39.851584,
                    "4d": -2.848386,
                    "4p": -5.769326,
                    "4s": -7.455966,
                    "5p": -0.504903,
                    "5s": -0.915819,
                    "6s": -0.078699
                },
                "Atomic radius": 2.6,
                "Atomic radius calculated": 2.98,
                "Boiling point": "944 K",
                "Brinell hardness": "0.14 MN m<sup>-2</sup>",
                "Bulk modulus": "1.6 GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    1
                ],
                "Critical temperature": "1938 K",
                "Density of solid": "1879 kg m<sup>-3</sup>",
                "Electrical resistivity": "21 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].6s<sup>1</sup>",
                "ICSD oxidation states": [
                    1
                ],
                "Ionic radii": {
                    "1": 1.81
                },
                "Liquid range": "642.41 K",
                "Melting point": "301.59 K",
                "Mendeleev no": 8,
                "Mineral hardness": "0.2",
                "Molar volume": "70.94 cm<sup>3</sup>",
                "Name": "Cesium",
                "Oxidation states": [
                    1
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "1": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.81,
                                "ionic_radius": 1.67
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.88,
                                "ionic_radius": 1.74
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.92,
                                "ionic_radius": 1.78
                            }
                        },
                        "X": {
                            "": {
                                "crystal_radius": 1.95,
                                "ionic_radius": 1.81
                            }
                        },
                        "XI": {
                            "": {
                                "crystal_radius": 1.99,
                                "ionic_radius": 1.85
                            }
                        },
                        "XII": {
                            "": {
                                "crystal_radius": 2.02,
                                "ionic_radius": 1.88
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "36 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 3.43,
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 0.79,
                "Youngs modulus": "1.7 GPa",
                "Metallic radius": 2.719,
                "iupac_ordering": 7,
                "IUPAC ordering": 7
            },
            "Cu": {
                "Atomic mass": 63.546,
                "Atomic no": 29,
                "Atomic orbitals": {
                    "1s": -320.78852,
                    "2p": -33.481247,
                    "2s": -38.14131,
                    "3d": -0.202272,
                    "3p": -2.609244,
                    "3s": -4.057453,
                    "4s": -0.172056
                },
                "Atomic radius": 1.35,
                "Atomic radius calculated": 1.45,
                "Boiling point": "3200 K",
                "Brinell hardness": "874 MN m<sup>-2</sup>",
                "Bulk modulus": "140 GPa",
                "Coefficient of linear thermal expansion": "16.5 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    2
                ],
                "Critical temperature": "no data K",
                "Density of solid": "8920 kg m<sup>-3</sup>",
                "Electrical resistivity": "1.72 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ar].3d<sup>10</sup>.4s<sup>1</sup>",
                "ICSD oxidation states": [
                    1,
                    2,
                    3
                ],
                "Ionic radii": {
                    "1": 0.91,
                    "2": 0.87,
                    "3": 0.68
                },
                "Liquid range": "1842.23 K",
                "Melting point": "1357.77 K",
                "Mendeleev no": 72,
                "Mineral hardness": "3.0",
                "Molar volume": "7.11 cm<sup>3</sup>",
                "Name": "Copper",
                "Oxidation states": [
                    1,
                    2,
                    3,
                    4
                ],
                "Poissons ratio": "0.34",
                "Reflectivity": "90 %",
                "Refractive index": "no data",
                "Rigidity modulus": "48 GPa",
                "Shannon radii": {
                    "1": {
                        "II": {
                            "": {
                                "crystal_radius": 0.6,
                                "ionic_radius": 0.46
                            }
                        },
                        "IV": {
                            "": {
                                "crystal_radius": 0.74,
                                "ionic_radius": 0.6
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.91,
                                "ionic_radius": 0.77
                            }
                        }
                    },
                    "2": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.71,
                                "ionic_radius": 0.57
                            }
                        },
                        "IVSQ": {
                            "": {
                                "crystal_radius": 0.71,
                                "ionic_radius": 0.57
                            }
                        },
                        "V": {
                            "": {
                                "crystal_radius": 0.79,
                                "ionic_radius": 0.65
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.87,
                                "ionic_radius": 0.73
                            }
                        }
                    },
                    "3": {
                        "VI": {
                            "Low Spin": {
                                "crystal_radius": 0.68,
                                "ionic_radius": 0.54
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "400 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.4,
                "Velocity of sound": "3570 m s<sup>-1</sup>",
                "Vickers hardness": "369 MN m<sup>-2</sup>",
                "X": 1.9,
                "Youngs modulus": "130 GPa",
                "NMR Quadrupole Moment": {
                    "Cu-63": -220.15,
                    "Cu-65": -204.14
                },
                "Metallic radius": 1.278,
                "iupac_ordering": 73,
                "IUPAC ordering": 73
            },
            "Dy": {
                "Atomic mass": 162.5,
                "Atomic no": 66,
                "Atomic orbitals": {
                    "1s": -1843.229585,
                    "2p": -281.558531,
                    "2s": -295.342856,
                    "3d": -47.4867,
                    "3p": -59.091931,
                    "3s": -65.299442,
                    "4d": -5.686352,
                    "4f": -0.265302,
                    "4p": -10.094091,
                    "4s": -12.551251,
                    "5p": -0.90349,
                    "5s": -1.547977,
                    "6s": -0.132769
                },
                "Atomic radius": 1.75,
                "Atomic radius calculated": 2.28,
                "Boiling point": "2840 K",
                "Brinell hardness": "500 MN m<sup>-2</sup>",
                "Bulk modulus": "41 GPa",
                "Coefficient of linear thermal expansion": "9.9 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "8551 kg m<sup>-3</sup>",
                "Electrical resistivity": "92.6 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>10</sup>.6s<sup>2</sup>",
                "ICSD oxidation states": [
                    3
                ],
                "Ionic radii": {
                    "2": 1.21,
                    "3": 1.052
                },
                "Liquid range": "1160 K",
                "Melting point": "1680 K",
                "Mendeleev no": 24,
                "Mineral hardness": "no data",
                "Molar volume": "19.01 cm<sup>3</sup>",
                "Name": "Dysprosium",
                "Oxidation states": [
                    2,
                    3
                ],
                "Poissons ratio": "0.25",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "25 GPa",
                "Shannon radii": {
                    "2": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.21,
                                "ionic_radius": 1.07
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.27,
                                "ionic_radius": 1.13
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.33,
                                "ionic_radius": 1.19
                            }
                        }
                    },
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.052,
                                "ionic_radius": 0.912
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.11,
                                "ionic_radius": 0.97
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.167,
                                "ionic_radius": 1.027
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.223,
                                "ionic_radius": 1.083
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "11 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "2710 m s<sup>-1</sup>",
                "Vickers hardness": "540 MN m<sup>-2</sup>",
                "X": 1.22,
                "Youngs modulus": "61 GPa",
                "Metallic radius": 1.773,
                "iupac_ordering": 38,
                "IUPAC ordering": 38
            },
            "Er": {
                "Atomic mass": 167.259,
                "Atomic no": 68,
                "Atomic orbitals": {
                    "1s": -1961.799176,
                    "2p": -302.01827,
                    "2s": -316.310631,
                    "3d": -51.682149,
                    "3p": -63.818655,
                    "3s": -70.310142,
                    "4d": -6.127443,
                    "4f": -0.278577,
                    "4p": -10.819574,
                    "4s": -13.423547,
                    "5p": -0.935202,
                    "5s": -1.616073,
                    "6s": -0.134905
                },
                "Atomic radius": 1.75,
                "Atomic radius calculated": 2.26,
                "Boiling point": "3141 K",
                "Brinell hardness": "814 MN m<sup>-2</sup>",
                "Bulk modulus": "44 GPa",
                "Coefficient of linear thermal expansion": "12.2 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "9066 kg m<sup>-3</sup>",
                "Electrical resistivity": "86.0 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>12</sup>.6s<sup>2</sup>",
                "ICSD oxidation states": [
                    3
                ],
                "Ionic radii": {
                    "3": 1.03
                },
                "Liquid range": "1371 K",
                "Melting point": "1802 K",
                "Mendeleev no": 22,
                "Mineral hardness": "no data",
                "Molar volume": "18.46 cm<sup>3</sup>",
                "Name": "Erbium",
                "Oxidation states": [
                    3
                ],
                "Poissons ratio": "0.24",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "28 GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.03,
                                "ionic_radius": 0.89
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.085,
                                "ionic_radius": 0.945
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.144,
                                "ionic_radius": 1.004
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.202,
                                "ionic_radius": 1.062
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "15 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "2830 m s<sup>-1</sup>",
                "Vickers hardness": "589 MN m<sup>-2</sup>",
                "X": 1.24,
                "Youngs modulus": "70 GPa",
                "Metallic radius": 1.756,
                "iupac_ordering": 36,
                "IUPAC ordering": 36
            },
            "Es": {
                "Atomic mass": 252.0,
                "Atomic no": 99,
                "Atomic orbitals": "no data",
                "Atomic radius": "no data",
                "Atomic radius calculated": "no data",
                "Boiling point": "no data K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "no data kg m<sup>-3</sup>",
                "Electrical resistivity": "no data 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Rn].5f<sup>11</sup>.7s<sup>2</sup>",
                "Liquid range": "no data K",
                "Melting point": "1133 K",
                "Mendeleev no": 38,
                "Mineral hardness": "no data",
                "Molar volume": "28.52 cm<sup>3</sup>",
                "Name": "Einsteinium",
                "Oxidation states": [
                    2,
                    3
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "no data W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.3,
                "Youngs modulus": "no data GPa",
                "Metallic radius": 1.86,
                "iupac_ordering": 22,
                "IUPAC ordering": 22
            },
            "Eu": {
                "Atomic mass": 151.964,
                "Atomic no": 63,
                "Atomic orbitals": {
                    "1s": -1672.309322,
                    "2p": -252.176697,
                    "2s": -265.199534,
                    "3d": -41.465518,
                    "3p": -52.281987,
                    "3s": -58.068128,
                    "4d": -5.03242,
                    "4f": -0.232773,
                    "4p": -9.025455,
                    "4s": -11.267747,
                    "5p": -0.853575,
                    "5s": -1.444087,
                    "6s": -0.129426
                },
                "Atomic radius": 1.85,
                "Atomic radius calculated": 2.31,
                "Boiling point": "1800 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "8.3 GPa",
                "Coefficient of linear thermal expansion": "35 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    2,
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "5244 kg m<sup>-3</sup>",
                "Electrical resistivity": "90 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>7</sup>.6s<sup>2</sup>",
                "ICSD oxidation states": [
                    2,
                    3
                ],
                "Ionic radii": {
                    "2": 1.31,
                    "3": 1.087
                },
                "Liquid range": "701 K",
                "Melting point": "1099 K",
                "Mendeleev no": 18,
                "Mineral hardness": "no data",
                "Molar volume": "28.97 cm<sup>3</sup>",
                "Name": "Europium",
                "Oxidation states": [
                    2,
                    3
                ],
                "Poissons ratio": "0.15",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "7.9 GPa",
                "Shannon radii": {
                    "2": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.31,
                                "ionic_radius": 1.17
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.34,
                                "ionic_radius": 1.2
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.39,
                                "ionic_radius": 1.25
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.44,
                                "ionic_radius": 1.3
                            }
                        },
                        "X": {
                            "": {
                                "crystal_radius": 1.49,
                                "ionic_radius": 1.35
                            }
                        }
                    },
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.087,
                                "ionic_radius": 0.947
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.15,
                                "ionic_radius": 1.01
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.206,
                                "ionic_radius": 1.066
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.26,
                                "ionic_radius": 1.12
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "14 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "167 MN m<sup>-2</sup>",
                "X": 1.2,
                "Youngs modulus": "18 GPa",
                "Metallic radius": 2.041,
                "iupac_ordering": 41,
                "IUPAC ordering": 41
            },
            "F": {
                "Atomic mass": 18.9984032,
                "Atomic no": 9,
                "Atomic orbitals": {
                    "1s": -24.189391,
                    "2p": -0.415606,
                    "2s": -1.086859
                },
                "Atomic radius": 0.5,
                "Atomic radius calculated": 0.42,
                "Boiling point": "85.03 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    -1
                ],
                "Critical temperature": "144 K",
                "Density of solid": "no data kg m<sup>-3</sup>",
                "Electrical resistivity": "no data 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[He].2s<sup>2</sup>.2p<sup>5</sup>",
                "ICSD oxidation states": [
                    -1
                ],
                "Ionic radii": {
                    "-1": 1.19,
                    "7": 0.22
                },
                "Liquid range": "31.5 K",
                "Melting point": "53.53 K",
                "Mendeleev no": 102,
                "Mineral hardness": "no data",
                "Molar volume": "11.20 cm<sup>3</sup>",
                "Name": "Fluorine",
                "Oxidation states": [
                    -1
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "1.000195",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "-1": {
                        "II": {
                            "": {
                                "crystal_radius": 1.145,
                                "ionic_radius": 1.285
                            }
                        },
                        "III": {
                            "": {
                                "crystal_radius": 1.16,
                                "ionic_radius": 1.3
                            }
                        },
                        "IV": {
                            "": {
                                "crystal_radius": 1.17,
                                "ionic_radius": 1.31
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 1.19,
                                "ionic_radius": 1.33
                            }
                        }
                    },
                    "7": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.22,
                                "ionic_radius": 0.08
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "0.0277 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.47,
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 3.98,
                "Youngs modulus": "no data GPa",
                "NMR Quadrupole Moment": {
                    "F-19": -94.2
                },
                "Metallic radius": "no data",
                "iupac_ordering": 102,
                "IUPAC ordering": 102
            },
            "Fe": {
                "Atomic mass": 55.845,
                "Atomic no": 26,
                "Atomic orbitals": {
                    "1s": -254.225505,
                    "2p": -25.551766,
                    "2s": -29.56486,
                    "3d": -0.295049,
                    "3p": -2.187523,
                    "3s": -3.360621,
                    "4s": -0.197978
                },
                "Atomic radius": 1.4,
                "Atomic radius calculated": 1.56,
                "Boiling point": "3134 K",
                "Brinell hardness": "490 MN m<sup>-2</sup>",
                "Bulk modulus": "170 GPa",
                "Coefficient of linear thermal expansion": "11.8 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    2,
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "7874 kg m<sup>-3</sup>",
                "Electrical resistivity": "10 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ar].3d<sup>6</sup>.4s<sup>2</sup>",
                "ICSD oxidation states": [
                    2,
                    3
                ],
                "Ionic radii": {
                    "2": 0.92,
                    "3": 0.785
                },
                "Ionic radii hs": {
                    "2": 0.92,
                    "3": 0.785
                },
                "Ionic radii ls": {
                    "2": 0.75,
                    "3": 0.69,
                    "4": 0.725,
                    "6": 0.39
                },
                "Liquid range": "1323 K",
                "Melting point": "1811 K",
                "Mendeleev no": 61,
                "Mineral hardness": "4.0",
                "Molar volume": "7.09 cm<sup>3</sup>",
                "Name": "Iron",
                "Oxidation states": [
                    -2,
                    -1,
                    1,
                    2,
                    3,
                    4,
                    5,
                    6
                ],
                "Poissons ratio": "0.29",
                "Reflectivity": "65 %",
                "Refractive index": "no data",
                "Rigidity modulus": "82 GPa",
                "Shannon radii": {
                    "2": {
                        "IV": {
                            "High Spin": {
                                "crystal_radius": 0.77,
                                "ionic_radius": 0.63
                            }
                        },
                        "IVSQ": {
                            "High Spin": {
                                "crystal_radius": 0.78,
                                "ionic_radius": 0.64
                            }
                        },
                        "VI": {
                            "Low Spin": {
                                "crystal_radius": 0.75,
                                "ionic_radius": 0.61
                            },
                            "High Spin": {
                                "crystal_radius": 0.92,
                                "ionic_radius": 0.78
                            }
                        },
                        "VIII": {
                            "High Spin": {
                                "crystal_radius": 1.06,
                                "ionic_radius": 0.92
                            }
                        }
                    },
                    "3": {
                        "IV": {
                            "High Spin": {
                                "crystal_radius": 0.63,
                                "ionic_radius": 0.49
                            }
                        },
                        "V": {
                            "": {
                                "crystal_radius": 0.72,
                                "ionic_radius": 0.58
                            }
                        },
                        "VI": {
                            "Low Spin": {
                                "crystal_radius": 0.69,
                                "ionic_radius": 0.55
                            },
                            "High Spin": {
                                "crystal_radius": 0.785,
                                "ionic_radius": 0.645
                            }
                        },
                        "VIII": {
                            "High Spin": {
                                "crystal_radius": 0.92,
                                "ionic_radius": 0.78
                            }
                        }
                    },
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.725,
                                "ionic_radius": 0.585
                            }
                        }
                    },
                    "6": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.39,
                                "ionic_radius": 0.25
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "80 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "4910 m s<sup>-1</sup>",
                "Vickers hardness": "608 MN m<sup>-2</sup>",
                "X": 1.83,
                "Youngs modulus": "211 GPa",
                "NMR Quadrupole Moment": {
                    "Fe-57": 160.0
                },
                "Metallic radius": 1.277,
                "iupac_ordering": 64,
                "IUPAC ordering": 64
            },
            "Fm": {
                "Atomic mass": 257.0,
                "Atomic no": 100,
                "Atomic orbitals": "no data",
                "Atomic radius": "no data",
                "Atomic radius calculated": "no data",
                "Boiling point": "no data K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "no data kg m<sup>-3</sup>",
                "Electrical resistivity": "no data 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Rn].5f<sup>12</sup>.7s<sup>2</sup>",
                "Liquid range": "no data K",
                "Melting point": "about 1800 K",
                "Mendeleev no": 37,
                "Mineral hardness": "no data",
                "Molar volume": "no data cm<sup>3</sup>",
                "Name": "Fermium",
                "Oxidation states": [
                    2,
                    3
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "no data W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.3,
                "Youngs modulus": "no data GPa",
                "Metallic radius": "no data",
                "iupac_ordering": 21,
                "IUPAC ordering": 21
            },
            "Fr": {
                "Atomic mass": 223.0,
                "Atomic no": 87,
                "Atomic orbitals": {
                    "1s": -3283.263399,
                    "2p": -542.41424,
                    "2s": -561.73045,
                    "3d": -111.085223,
                    "3p": -128.607136,
                    "3s": -137.959632,
                    "4d": -20.812462,
                    "4f": -10.050648,
                    "4p": -28.648131,
                    "4s": -32.861013,
                    "5d": -2.360991,
                    "5p": -4.97328,
                    "5s": -6.509516,
                    "6p": -0.466197,
                    "6s": -0.841848,
                    "7s": -0.076176
                },
                "Atomic radius": "no data",
                "Atomic radius calculated": "no data",
                "Boiling point": "no data K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    1
                ],
                "Critical temperature": "no data K",
                "Density of solid": "no data kg m<sup>-3</sup>",
                "Electrical resistivity": "no data 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Rn].7s<sup>1</sup>",
                "Ionic radii": {
                    "1": 1.94
                },
                "Liquid range": "no data K",
                "Melting point": "maybe about 300 K",
                "Mendeleev no": 7,
                "Mineral hardness": "no data",
                "Molar volume": "no data cm<sup>3</sup>",
                "Name": "Francium",
                "Oxidation states": [
                    1
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "1": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.94,
                                "ionic_radius": 1.8
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "no data W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 3.48,
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 0.7,
                "Youngs modulus": "no data GPa",
                "Metallic radius": "no data",
                "iupac_ordering": 6,
                "IUPAC ordering": 6
            },
            "Ga": {
                "Atomic mass": 69.723,
                "Atomic no": 31,
                "Atomic orbitals": {
                    "1s": -370.170639,
                    "2p": -40.093339,
                    "2s": -45.200869,
                    "3d": -0.736204,
                    "3p": -3.584666,
                    "3s": -5.241645,
                    "4p": -0.101634,
                    "4s": -0.328019
                },
                "Atomic radius": 1.3,
                "Atomic radius calculated": 1.36,
                "Boiling point": "2477 K",
                "Brinell hardness": "60 MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "120 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "5904 kg m<sup>-3</sup>",
                "Electrical resistivity": "about  14 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ar].3d<sup>10</sup>.4s<sup>2</sup>.4p<sup>1</sup>",
                "ICSD oxidation states": [
                    2,
                    3
                ],
                "Ionic radii": {
                    "3": 0.76
                },
                "Liquid range": "2174.09 K",
                "Melting point": "302.91 K",
                "Mendeleev no": 81,
                "Mineral hardness": "1.5",
                "Molar volume": "11.80 cm<sup>3</sup>",
                "Name": "Gallium",
                "Oxidation states": [
                    1,
                    2,
                    3
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "3": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.61,
                                "ionic_radius": 0.47
                            }
                        },
                        "V": {
                            "": {
                                "crystal_radius": 0.69,
                                "ionic_radius": 0.55
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.76,
                                "ionic_radius": 0.62
                            }
                        }
                    }
                },
                "Superconduction temperature": "1.083 K",
                "Thermal conductivity": "29 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.87,
                "Velocity of sound": "2740 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.81,
                "Youngs modulus": "no data GPa",
                "Metallic radius": 1.35,
                "iupac_ordering": 79,
                "IUPAC ordering": 79
            },
            "Gd": {
                "Atomic mass": 157.25,
                "Atomic no": 64,
                "Atomic orbitals": {
                    "1s": -1728.625195,
                    "2p": -262.081616,
                    "2s": -275.36313,
                    "3d": -43.754556,
                    "3p": -54.836922,
                    "3s": -60.764408,
                    "4d": -5.531835,
                    "4f": -0.489012,
                    "4p": -9.669866,
                    "4s": -11.986486,
                    "5d": -0.12722,
                    "5p": -0.978749,
                    "5s": -1.608477,
                    "6s": -0.143627
                },
                "Atomic radius": 1.8,
                "Atomic radius calculated": 2.33,
                "Boiling point": "3523 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "38 GPa",
                "Coefficient of linear thermal expansion": "9.4 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "7901 kg m<sup>-3</sup>",
                "Electrical resistivity": "131 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>7</sup>.5d<sup>1</sup>.6s<sup>2</sup>",
                "ICSD oxidation states": [
                    3
                ],
                "Ionic radii": {
                    "3": 1.075
                },
                "Liquid range": "1938 K",
                "Melting point": "1585 K",
                "Mendeleev no": 27,
                "Mineral hardness": "no data",
                "Molar volume": "19.90 cm<sup>3</sup>",
                "Name": "Gadolinium",
                "Oxidation states": [
                    1,
                    2,
                    3
                ],
                "Poissons ratio": "0.26",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "22 GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.078,
                                "ionic_radius": 0.938
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.14,
                                "ionic_radius": 1.0
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.193,
                                "ionic_radius": 1.053
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.247,
                                "ionic_radius": 1.107
                            }
                        }
                    }
                },
                "Superconduction temperature": "1.083 K",
                "Thermal conductivity": "11 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "2680 m s<sup>-1</sup>",
                "Vickers hardness": "570 MN m<sup>-2</sup>",
                "X": 1.2,
                "Youngs modulus": "55 GPa",
                "Metallic radius": 1.802,
                "iupac_ordering": 40,
                "IUPAC ordering": 40
            },
            "Ge": {
                "Atomic mass": 72.64,
                "Atomic no": 32,
                "Atomic orbitals": {
                    "1s": -396.292991,
                    "2p": -43.720129,
                    "2s": -49.055282,
                    "3d": -1.117316,
                    "3p": -4.194822,
                    "3s": -5.961472,
                    "4p": -0.149882,
                    "4s": -0.426523
                },
                "Atomic radius": 1.25,
                "Atomic radius calculated": 1.25,
                "Boiling point": "3093 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "6 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    -4,
                    2,
                    4
                ],
                "Critical temperature": "no data K",
                "Density of solid": "5323 kg m<sup>-3</sup>",
                "Electrical resistivity": "about 50000 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ar].3d<sup>10</sup>.4s<sup>2</sup>.4p<sup>2</sup>",
                "ICSD oxidation states": [
                    2,
                    3,
                    4
                ],
                "Ionic radii": {
                    "2": 0.87,
                    "4": 0.67
                },
                "Liquid range": "1881.6 K",
                "Melting point": "1211.4 K",
                "Mendeleev no": 84,
                "Mineral hardness": "6.0",
                "Molar volume": "13.63 cm<sup>3</sup>",
                "Name": "Germanium",
                "Oxidation states": [
                    -4,
                    1,
                    2,
                    3,
                    4
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "2": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.87,
                                "ionic_radius": 0.73
                            }
                        }
                    },
                    "4": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.53,
                                "ionic_radius": 0.39
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.67,
                                "ionic_radius": 0.53
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "60 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 2.11,
                "Velocity of sound": "5400 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 2.01,
                "Youngs modulus": "no data GPa",
                "Metallic radius": 1.39,
                "iupac_ordering": 84,
                "IUPAC ordering": 84
            },
            "H": {
                "Atomic mass": 1.00794,
                "Atomic no": 1,
                "Atomic orbitals": {
                    "1s": -0.233471
                },
                "Atomic radius": 0.25,
                "Atomic radius calculated": 0.53,
                "Boiling point": "20.28 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    -1,
                    1
                ],
                "Critical temperature": "33 K",
                "Density of solid": "no data kg m<sup>-3</sup>",
                "Electrical resistivity": "no data 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "1s<sup>1</sup>",
                "ICSD oxidation states": [
                    1,
                    -1
                ],
                "Liquid range": "6.27 K",
                "Melting point": "14.01 K",
                "Mendeleev no": 103,
                "Mineral hardness": "no data",
                "Molar volume": "11.42 cm<sup>3</sup>",
                "Name": "Hydrogen",
                "Oxidation states": [
                    -1,
                    1
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "1.000132 (gas; liquid 1.12)(no units)",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "1": {
                        "I": {
                            "": {
                                "crystal_radius": -0.24,
                                "ionic_radius": -0.38
                            }
                        },
                        "II": {
                            "": {
                                "crystal_radius": -0.04,
                                "ionic_radius": -0.18
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "0.1805 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.2,
                "Velocity of sound": "1270 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 2.2,
                "Youngs modulus": "no data GPa",
                "NMR Quadrupole Moment": {
                    "H-2": 2.86
                },
                "Metallic radius": "no data",
                "iupac_ordering": 92,
                "IUPAC ordering": 92
            },
            "He": {
                "Atomic mass": 4.002602,
                "Atomic no": 2,
                "Atomic orbitals": {
                    "1s": -0.570425
                },
                "Atomic radius": "no data",
                "Atomic radius calculated": 0.31,
                "Boiling point": "4.22 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Critical temperature": "5.19 K",
                "Density of solid": "no data kg m<sup>-3</sup>",
                "Electrical resistivity": "no data 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "1s<sup>2</sup>",
                "Liquid range": "3.27 K",
                "Max oxidation state": 0.0,
                "Melting point": "0.95 K",
                "Mendeleev no": 1,
                "Min oxidation state": 0.0,
                "Mineral hardness": "no data",
                "Molar volume": "21.0 cm<sup>3</sup>",
                "Name": "Helium",
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "1.000035 (gas; liquid 1.028)(no units)",
                "Rigidity modulus": "no data GPa",
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "0.1513 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.4,
                "Velocity of sound": "970 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "Youngs modulus": "no data GPa",
                "Metallic radius": "no data",
                "iupac_ordering": 5,
                "IUPAC ordering": 5
            },
            "Hf": {
                "Atomic mass": 178.49,
                "Atomic no": 72,
                "Atomic orbitals": {
                    "1s": -2210.65199,
                    "2p": -345.687023,
                    "2s": -361.006527,
                    "3d": -61.231443,
                    "3p": -74.452656,
                    "3s": -81.522812,
                    "4d": -7.676638,
                    "4f": -0.871574,
                    "4p": -12.971211,
                    "4s": -15.883625,
                    "5d": -0.143805,
                    "5p": -1.246441,
                    "5s": -2.049828,
                    "6s": -0.166465
                },
                "Atomic radius": 1.55,
                "Atomic radius calculated": 2.08,
                "Boiling point": "4876 K",
                "Brinell hardness": "1700 MN m<sup>-2</sup>",
                "Bulk modulus": "110 GPa",
                "Coefficient of linear thermal expansion": "5.9 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    4
                ],
                "Critical temperature": "no data K",
                "Density of solid": "13310 kg m<sup>-3</sup>",
                "Electrical resistivity": "34 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>14</sup>.5d<sup>2</sup>.6s<sup>2</sup>",
                "ICSD oxidation states": [
                    4
                ],
                "Ionic radii": {
                    "4": 0.85
                },
                "Liquid range": "2370 K",
                "Melting point": "2506 K",
                "Mendeleev no": 50,
                "Mineral hardness": "5.5",
                "Molar volume": "13.44 cm<sup>3</sup>",
                "Name": "Hafnium",
                "Oxidation states": [
                    2,
                    3,
                    4
                ],
                "Poissons ratio": "0.37",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "30 GPa",
                "Shannon radii": {
                    "4": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.72,
                                "ionic_radius": 0.58
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.85,
                                "ionic_radius": 0.71
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 0.9,
                                "ionic_radius": 0.76
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 0.97,
                                "ionic_radius": 0.83
                            }
                        }
                    }
                },
                "Superconduction temperature": "0.128 K",
                "Thermal conductivity": "23 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "3010 m s<sup>-1</sup>",
                "Vickers hardness": "1760 MN m<sup>-2</sup>",
                "X": 1.3,
                "Youngs modulus": "78 GPa",
                "Metallic radius": 1.58,
                "iupac_ordering": 50,
                "IUPAC ordering": 50
            },
            "Hg": {
                "Atomic mass": 200.59,
                "Atomic no": 80,
                "Atomic orbitals": {
                    "1s": -2755.022637,
                    "2p": -443.848676,
                    "2s": -461.27864,
                    "3d": -84.845492,
                    "3p": -100.328031,
                    "3s": -108.597921,
                    "4d": -13.019221,
                    "4f": -4.110291,
                    "4p": -19.636187,
                    "4s": -23.222921,
                    "5d": -0.452552,
                    "5p": -2.261975,
                    "5s": -3.423486,
                    "6s": -0.205137
                },
                "Atomic radius": 1.5,
                "Atomic radius calculated": 1.71,
                "Boiling point": "629.88 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "25 GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    1,
                    2
                ],
                "Critical temperature": "1750 K",
                "Density of solid": "no data kg m<sup>-3</sup>",
                "Electrical resistivity": "96 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>14</sup>.5d<sup>10</sup>.6s<sup>2</sup>",
                "ICSD oxidation states": [
                    1,
                    2
                ],
                "Ionic radii": {
                    "1": 1.33,
                    "2": 1.16
                },
                "Liquid range": "395.56 K",
                "Melting point": "234.32 K",
                "Mendeleev no": 74,
                "Mineral hardness": "1.5",
                "Molar volume": "14.09 cm<sup>3</sup>",
                "Name": "Mercury",
                "Oxidation states": [
                    1,
                    2,
                    4
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "73 %",
                "Refractive index": "1.000933",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "1": {
                        "III": {
                            "": {
                                "crystal_radius": 1.11,
                                "ionic_radius": 0.97
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 1.33,
                                "ionic_radius": 1.19
                            }
                        }
                    },
                    "2": {
                        "II": {
                            "": {
                                "crystal_radius": 0.83,
                                "ionic_radius": 0.69
                            }
                        },
                        "IV": {
                            "": {
                                "crystal_radius": 1.1,
                                "ionic_radius": 0.96
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 1.16,
                                "ionic_radius": 1.02
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.28,
                                "ionic_radius": 1.14
                            }
                        }
                    }
                },
                "Superconduction temperature": "3.95 K",
                "Thermal conductivity": "8.3 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.55,
                "Velocity of sound": "1407 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 2.0,
                "Youngs modulus": "no data GPa",
                "NMR Quadrupole Moment": {
                    "Hg-201": 387.6
                },
                "Metallic radius": 1.51,
                "iupac_ordering": 74,
                "IUPAC ordering": 74
            },
            "Ho": {
                "Atomic mass": 164.93032,
                "Atomic no": 67,
                "Atomic orbitals": {
                    "1s": -1902.051908,
                    "2p": -291.700994,
                    "2s": -305.739294,
                    "3d": -49.565996,
                    "3p": -61.436304,
                    "3s": -67.785492,
                    "4d": -5.906195,
                    "4f": -0.272677,
                    "4p": -10.455303,
                    "4s": -12.985498,
                    "5p": -0.919463,
                    "5s": -1.582088,
                    "6s": -0.133845
                },
                "Atomic radius": 1.75,
                "Atomic radius calculated": "no data",
                "Boiling point": "2993 K",
                "Brinell hardness": "746 MN m<sup>-2</sup>",
                "Bulk modulus": "40 GPa",
                "Coefficient of linear thermal expansion": "11.2 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "8795 kg m<sup>-3</sup>",
                "Electrical resistivity": "81.4 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>11</sup>.6s<sup>2</sup>",
                "ICSD oxidation states": [
                    3
                ],
                "Ionic radii": {
                    "3": 1.041
                },
                "Liquid range": "1259 K",
                "Melting point": "1734 K",
                "Mendeleev no": 23,
                "Mineral hardness": "no data",
                "Molar volume": "18.74 cm<sup>3</sup>",
                "Name": "Holmium",
                "Oxidation states": [
                    3
                ],
                "Poissons ratio": "0.23",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "26 GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.041,
                                "ionic_radius": 0.901
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.155,
                                "ionic_radius": 1.015
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.212,
                                "ionic_radius": 1.072
                            }
                        },
                        "X": {
                            "": {
                                "crystal_radius": 1.26,
                                "ionic_radius": 1.12
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "16 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "2760 m s<sup>-1</sup>",
                "Vickers hardness": "481 MN m<sup>-2</sup>",
                "X": 1.23,
                "Youngs modulus": "65 GPa",
                "Metallic radius": 1.765,
                "iupac_ordering": 37,
                "IUPAC ordering": 37
            },
            "I": {
                "Atomic mass": 126.90447,
                "Atomic no": 53,
                "Atomic orbitals": {
                    "1s": -1161.787047,
                    "2p": -164.603788,
                    "2s": -175.073804,
                    "3d": -22.600693,
                    "3p": -30.831092,
                    "3s": -35.243351,
                    "4d": -1.938179,
                    "4p": -4.572522,
                    "4s": -6.115811,
                    "5p": -0.267904,
                    "5s": -0.596339
                },
                "Atomic radius": 1.4,
                "Atomic radius calculated": 1.15,
                "Boiling point": "457.4 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "7.7 GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    -1,
                    1,
                    3,
                    5,
                    7
                ],
                "Critical temperature": "819 K",
                "Density of solid": "4940 kg m<sup>-3</sup>",
                "Electrical resistivity": "&gt; 10<sup>15</sup>10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Kr].4d<sup>10</sup>.5s<sup>2</sup>.5p<sup>5</sup>",
                "ICSD oxidation states": [
                    5,
                    -1
                ],
                "Ionic radii": {
                    "-1": 2.06,
                    "5": 1.09,
                    "7": 0.67
                },
                "Liquid range": "70.55 K",
                "Melting point": "386.85 K",
                "Mendeleev no": 97,
                "Mineral hardness": "no data",
                "Molar volume": "25.72 cm<sup>3</sup>",
                "Name": "Iodine",
                "Oxidation states": [
                    -1,
                    1,
                    3,
                    5,
                    7
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "-1": {
                        "VI": {
                            "": {
                                "crystal_radius": 2.06,
                                "ionic_radius": 2.2
                            }
                        }
                    },
                    "5": {
                        "IIIPY": {
                            "": {
                                "crystal_radius": 0.58,
                                "ionic_radius": 0.44
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 1.09,
                                "ionic_radius": 0.95
                            }
                        }
                    },
                    "7": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.56,
                                "ionic_radius": 0.42
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.67,
                                "ionic_radius": 0.53
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "0.449 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.98,
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 2.66,
                "Youngs modulus": "no data GPa",
                "NMR Quadrupole Moment": {
                    "I-127": -696.12,
                    "I-129": -604.1
                },
                "Metallic radius": "no data",
                "iupac_ordering": 99,
                "IUPAC ordering": 99
            },
            "In": {
                "Atomic mass": 114.818,
                "Atomic no": 49,
                "Atomic orbitals": {
                    "1s": -983.647445,
                    "2p": -134.628845,
                    "2s": -144.078357,
                    "3d": -16.139823,
                    "3p": -23.345778,
                    "3s": -27.2206,
                    "4d": -0.730481,
                    "4p": -2.795832,
                    "4s": -4.062639,
                    "5p": -0.101782,
                    "5s": -0.290497
                },
                "Atomic radius": 1.55,
                "Atomic radius calculated": 1.56,
                "Boiling point": "2345 K",
                "Brinell hardness": "8.83 MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "32.1 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "7310 kg m<sup>-3</sup>",
                "Electrical resistivity": "8 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Kr].4d<sup>10</sup>.5s<sup>2</sup>.5p<sup>1</sup>",
                "ICSD oxidation states": [
                    1,
                    2,
                    3
                ],
                "Ionic radii": {
                    "3": 0.94
                },
                "Liquid range": "1915.25 K",
                "Melting point": "429.75 K",
                "Mendeleev no": 79,
                "Mineral hardness": "1.2",
                "Molar volume": "15.76 cm<sup>3</sup>",
                "Name": "Indium",
                "Oxidation states": [
                    1,
                    2,
                    3
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "3": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.76,
                                "ionic_radius": 0.62
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.94,
                                "ionic_radius": 0.8
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.06,
                                "ionic_radius": 0.92
                            }
                        }
                    }
                },
                "Superconduction temperature": "3.41 K",
                "Thermal conductivity": "82 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.93,
                "Velocity of sound": "1215 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.78,
                "Youngs modulus": "11 GPa",
                "NMR Quadrupole Moment": {
                    "In-113": 759.8,
                    "In-115": 770.8
                },
                "Metallic radius": 1.67,
                "iupac_ordering": 78,
                "IUPAC ordering": 78
            },
            "Ir": {
                "Atomic mass": 192.217,
                "Atomic no": 77,
                "Atomic orbitals": {
                    "1s": -2543.761342,
                    "2p": -405.526834,
                    "2s": -422.159424,
                    "3d": -75.485027,
                    "3p": -90.108427,
                    "3s": -97.923081,
                    "4d": -10.856593,
                    "4f": -2.738339,
                    "4p": -16.966578,
                    "4s": -20.29429,
                    "5d": -0.335189,
                    "5p": -1.883349,
                    "5s": -2.909174,
                    "6s": -0.195511
                },
                "Atomic radius": 1.35,
                "Atomic radius calculated": 1.8,
                "Boiling point": "4701 K",
                "Brinell hardness": "1670 MN m<sup>-2</sup>",
                "Bulk modulus": "320 GPa",
                "Coefficient of linear thermal expansion": "6.4 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3,
                    4
                ],
                "Critical temperature": "no data K",
                "Density of solid": "22650 kg m<sup>-3</sup>",
                "Electrical resistivity": "4.7 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>14</sup>.5d<sup>7</sup>.6s<sup>2</sup>",
                "ICSD oxidation states": [
                    3,
                    4,
                    5
                ],
                "Ionic radii": {
                    "3": 0.82,
                    "4": 0.765,
                    "5": 0.71
                },
                "Liquid range": "1962 K",
                "Melting point": "2739 K",
                "Mendeleev no": 66,
                "Mineral hardness": "6.5",
                "Molar volume": "8.52 cm<sup>3</sup>",
                "Name": "Iridium",
                "Oxidation states": [
                    -3,
                    -1,
                    1,
                    2,
                    3,
                    4,
                    5,
                    6
                ],
                "Poissons ratio": "0.26",
                "Reflectivity": "78 %",
                "Refractive index": "no data",
                "Rigidity modulus": "210 GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.82,
                                "ionic_radius": 0.68
                            }
                        }
                    },
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.765,
                                "ionic_radius": 0.625
                            }
                        }
                    },
                    "5": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.71,
                                "ionic_radius": 0.57
                            }
                        }
                    }
                },
                "Superconduction temperature": "0.11 K",
                "Thermal conductivity": "150 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "4825 m s<sup>-1</sup>",
                "Vickers hardness": "1760 MN m<sup>-2</sup>",
                "X": 2.2,
                "Youngs modulus": "528 GPa",
                "Metallic radius": 1.357,
                "iupac_ordering": 65,
                "IUPAC ordering": 65
            },
            "K": {
                "Atomic mass": 39.0983,
                "Atomic no": 19,
                "Atomic orbitals": {
                    "1s": -128.414957,
                    "2p": -10.283851,
                    "2s": -12.839001,
                    "3p": -0.693776,
                    "3s": -1.281897,
                    "4s": -0.088815
                },
                "Atomic radius": 2.2,
                "Atomic radius calculated": 2.43,
                "Boiling point": "1032 K",
                "Brinell hardness": "0.363 MN m<sup>-2</sup>",
                "Bulk modulus": "3.1 GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    1
                ],
                "Critical temperature": "2223 K",
                "Density of solid": "856 kg m<sup>-3</sup>",
                "Electrical resistivity": "7.5 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ar].4s<sup>1</sup>",
                "ICSD oxidation states": [
                    1
                ],
                "Ionic radii": {
                    "1": 1.52
                },
                "Liquid range": "695.47 K",
                "Melting point": "336.53 K",
                "Mendeleev no": 10,
                "Mineral hardness": "0.4",
                "Molar volume": "45.94 cm<sup>3</sup>",
                "Name": "Potassium",
                "Oxidation states": [
                    -1,
                    1
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "1.3 GPa",
                "Shannon radii": {
                    "1": {
                        "IV": {
                            "": {
                                "crystal_radius": 1.51,
                                "ionic_radius": 1.37
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 1.52,
                                "ionic_radius": 1.38
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.6,
                                "ionic_radius": 1.46
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.65,
                                "ionic_radius": 1.51
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.69,
                                "ionic_radius": 1.55
                            }
                        },
                        "X": {
                            "": {
                                "crystal_radius": 1.73,
                                "ionic_radius": 1.59
                            }
                        },
                        "XII": {
                            "": {
                                "crystal_radius": 1.78,
                                "ionic_radius": 1.64
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "100 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 2.75,
                "Velocity of sound": "2000 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 0.82,
                "Youngs modulus": "no data GPa",
                "NMR Quadrupole Moment": {
                    "K-39": 58.5,
                    "K-40": -73.0,
                    "K-41": 71.1
                },
                "Metallic radius": 2.381,
                "iupac_ordering": 9,
                "IUPAC ordering": 9
            },
            "Kr": {
                "Atomic mass": 83.798,
                "Atomic no": 36,
                "Atomic orbitals": {
                    "1s": -509.982989,
                    "2p": -60.017328,
                    "2s": -66.285953,
                    "3d": -3.074109,
                    "3p": -7.086634,
                    "3s": -9.315192,
                    "4p": -0.34634,
                    "4s": -0.820574
                },
                "Atomic radius": "no data",
                "Atomic radius calculated": 0.88,
                "Boiling point": "119.93 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Critical temperature": "209.4 K",
                "Density of solid": "no data kg m<sup>-3</sup>",
                "Electrical resistivity": "no data 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ar].3d<sup>10</sup>.4s<sup>2</sup>.4p<sup>6</sup>",
                "Liquid range": "4.14 K",
                "Max oxidation state": 0.0,
                "Melting point": "115.79 K",
                "Mendeleev no": 4,
                "Min oxidation state": 0.0,
                "Mineral hardness": "no data",
                "Molar volume": "27.99 cm<sup>3</sup>",
                "Name": "Krypton",
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "1.000427",
                "Rigidity modulus": "no data GPa",
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "0.00943 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 2.02,
                "Velocity of sound": "1120 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 3.0,
                "Youngs modulus": "no data GPa",
                "Metallic radius": "no data",
                "iupac_ordering": 2,
                "IUPAC ordering": 2
            },
            "La": {
                "Atomic mass": 138.90547,
                "Atomic no": 57,
                "Atomic orbitals": {
                    "1s": -1355.622446,
                    "2p": -198.325243,
                    "2s": -209.831151,
                    "3d": -30.626696,
                    "3p": -39.895838,
                    "3s": -44.856283,
                    "4d": -3.95801,
                    "4p": -7.167724,
                    "4s": -9.000543,
                    "5d": -0.141085,
                    "5p": -0.824498,
                    "5s": -1.324936,
                    "6s": -0.132233
                },
                "Atomic radius": 1.95,
                "Atomic radius calculated": "no data",
                "Boiling point": "3743 K",
                "Brinell hardness": "363 MN m<sup>-2</sup>",
                "Bulk modulus": "28 GPa",
                "Coefficient of linear thermal expansion": "12.1 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "6146 kg m<sup>-3</sup>",
                "Electrical resistivity": "61.5 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].5d<sup>1</sup>.6s<sup>2</sup>",
                "ICSD oxidation states": [
                    2,
                    3
                ],
                "Ionic radii": {
                    "3": 1.172
                },
                "Liquid range": "2550 K",
                "Melting point": "1193 K",
                "Mendeleev no": 33,
                "Mineral hardness": "2.5",
                "Molar volume": "22.39 cm<sup>3</sup>",
                "Name": "Lanthanum",
                "Oxidation states": [
                    2,
                    3
                ],
                "Poissons ratio": "0.28",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "14 GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.172,
                                "ionic_radius": 1.032
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.24,
                                "ionic_radius": 1.1
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.3,
                                "ionic_radius": 1.16
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.356,
                                "ionic_radius": 1.216
                            }
                        },
                        "X": {
                            "": {
                                "crystal_radius": 1.41,
                                "ionic_radius": 1.27
                            }
                        },
                        "XII": {
                            "": {
                                "crystal_radius": 1.5,
                                "ionic_radius": 1.36
                            }
                        }
                    }
                },
                "Superconduction temperature": "6.00 K",
                "Thermal conductivity": "13 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "2475 m s<sup>-1</sup>",
                "Vickers hardness": "491 MN m<sup>-2</sup>",
                "X": 1.1,
                "Youngs modulus": "37 GPa",
                "NMR Quadrupole Moment": {
                    "La-139": 200.6
                },
                "Metallic radius": 1.877,
                "iupac_ordering": 47,
                "IUPAC ordering": 47
            },
            "Li": {
                "Atomic mass": 6.941,
                "Atomic no": 3,
                "Atomic orbitals": {
                    "1s": -1.878564,
                    "2s": -0.10554
                },
                "Atomic radius": 1.45,
                "Atomic radius calculated": 1.67,
                "Boiling point": "1615 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "11 GPa",
                "Coefficient of linear thermal expansion": "46 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    1
                ],
                "Critical temperature": "3223 K",
                "Density of solid": "535 kg m<sup>-3</sup>",
                "Electrical resistivity": "9.5 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[He].2s<sup>1</sup>",
                "ICSD oxidation states": [
                    1
                ],
                "Ionic radii": {
                    "1": 0.9
                },
                "Liquid range": "1161.31 K",
                "Melting point": "453.69 K",
                "Mendeleev no": 12,
                "Mineral hardness": "0.6",
                "Molar volume": "13.02 cm<sup>3</sup>",
                "Name": "Lithium",
                "Oxidation states": [
                    1
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "4.2 GPa",
                "Shannon radii": {
                    "1": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.73,
                                "ionic_radius": 0.59
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.9,
                                "ionic_radius": 0.76
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.06,
                                "ionic_radius": 0.92
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "85 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.82,
                "Velocity of sound": "6000 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 0.98,
                "Youngs modulus": "4.9 GPa",
                "NMR Quadrupole Moment": {
                    "Li-6": -0.808,
                    "Li-7": -40.1
                },
                "Metallic radius": 1.52,
                "iupac_ordering": 11,
                "IUPAC ordering": 11
            },
            "Lr": {
                "Atomic mass": 262.0,
                "Atomic no": 103,
                "Atomic orbitals": "no data",
                "Atomic radius": "no data",
                "Atomic radius calculated": "no data",
                "Boiling point": "no data K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "no data kg m<sup>-3</sup>",
                "Electrical resistivity": "no data 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Rn].5f<sup>14</sup>.7s<sup>2</sup>.7p<sup>1</sup> (tentative)",
                "Liquid range": "no data K",
                "Melting point": "about 1900 K",
                "Mendeleev no": 34,
                "Mineral hardness": "no data",
                "Molar volume": "no data cm<sup>3</sup>",
                "Name": "Lawrencium",
                "Oxidation states": [
                    3
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "no data W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.3,
                "Youngs modulus": "no data GPa",
                "Metallic radius": "no data",
                "iupac_ordering": 18,
                "IUPAC ordering": 18
            },
            "Lu": {
                "Atomic mass": 174.967,
                "Atomic no": 71,
                "Atomic orbitals": {
                    "1s": -2146.885351,
                    "2p": -334.330902,
                    "2s": -349.390492,
                    "3d": -58.592982,
                    "3p": -71.538779,
                    "3s": -78.462398,
                    "4d": -7.113364,
                    "4f": -0.568096,
                    "4p": -12.250904,
                    "4s": -15.08337,
                    "5d": -0.103686,
                    "5p": -1.111991,
                    "5s": -1.872086,
                    "6s": -0.155112
                },
                "Atomic radius": 1.75,
                "Atomic radius calculated": 2.17,
                "Boiling point": "3675 K",
                "Brinell hardness": "893 MN m<sup>-2</sup>",
                "Bulk modulus": "48 GPa",
                "Coefficient of linear thermal expansion": "9.9 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "9841 kg m<sup>-3</sup>",
                "Electrical resistivity": "58 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>14</sup>.5d<sup>1</sup>.6s<sup>2</sup>",
                "ICSD oxidation states": [
                    3
                ],
                "Ionic radii": {
                    "3": 1.001
                },
                "Liquid range": "1750 K",
                "Melting point": "1925 K",
                "Mendeleev no": 20,
                "Mineral hardness": "no data",
                "Molar volume": "17.78 cm<sup>3</sup>",
                "Name": "Lutetium",
                "Oxidation states": [
                    3
                ],
                "Poissons ratio": "0.26",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "27 GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.001,
                                "ionic_radius": 0.861
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.117,
                                "ionic_radius": 0.977
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.172,
                                "ionic_radius": 1.032
                            }
                        }
                    }
                },
                "Superconduction temperature": "0.022 K",
                "Thermal conductivity": "16 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "1160 MN m<sup>-2</sup>",
                "X": 1.27,
                "Youngs modulus": "69 GPa",
                "Metallic radius": 1.735,
                "iupac_ordering": 33,
                "IUPAC ordering": 33
            },
            "Md": {
                "Atomic mass": 258.0,
                "Atomic no": 101,
                "Atomic orbitals": "no data",
                "Atomic radius": "no data",
                "Atomic radius calculated": "no data",
                "Boiling point": "no data K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "no data kg m<sup>-3</sup>",
                "Electrical resistivity": "no data 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Rn].5f<sup>13</sup>.7s<sup>2</sup>",
                "Liquid range": "no data K",
                "Melting point": "about 1100 K",
                "Mendeleev no": 36,
                "Mineral hardness": "no data",
                "Molar volume": "no data cm<sup>3</sup>",
                "Name": "Mendelevium",
                "Oxidation states": [
                    2,
                    3
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "no data W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.3,
                "Youngs modulus": "no data GPa",
                "Metallic radius": "no data",
                "iupac_ordering": 20,
                "IUPAC ordering": 20
            },
            "Mg": {
                "Atomic mass": 24.305,
                "Atomic no": 12,
                "Atomic orbitals": {
                    "1s": -45.973167,
                    "2p": -1.71897,
                    "2s": -2.903746,
                    "3s": -0.175427
                },
                "Atomic radius": 1.5,
                "Atomic radius calculated": 1.45,
                "Boiling point": "1363 K",
                "Brinell hardness": "260 MN m<sup>-2</sup>",
                "Bulk modulus": "45 GPa",
                "Coefficient of linear thermal expansion": "8.2 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    2
                ],
                "Critical temperature": "no data K",
                "Density of solid": "1738 kg m<sup>-3</sup>",
                "Electrical resistivity": "4.4 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ne].3s<sup>2</sup>",
                "ICSD oxidation states": [
                    2
                ],
                "Ionic radii": {
                    "2": 0.86
                },
                "Liquid range": "440 K",
                "Melting point": "923 K",
                "Mendeleev no": 73,
                "Mineral hardness": "2.5",
                "Molar volume": "14.00 cm<sup>3</sup>",
                "Name": "Magnesium",
                "Oxidation states": [
                    1,
                    2
                ],
                "Poissons ratio": "0.29",
                "Reflectivity": "74 %",
                "Refractive index": "no data",
                "Rigidity modulus": "17 GPa",
                "Shannon radii": {
                    "2": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.71,
                                "ionic_radius": 0.57
                            }
                        },
                        "V": {
                            "": {
                                "crystal_radius": 0.8,
                                "ionic_radius": 0.66
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.86,
                                "ionic_radius": 0.72
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.03,
                                "ionic_radius": 0.89
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "160 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.73,
                "Velocity of sound": "4602 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.31,
                "Youngs modulus": "45 GPa",
                "NMR Quadrupole Moment": {
                    "Mg-25": 199.4
                },
                "Metallic radius": 1.6,
                "iupac_ordering": 16,
                "IUPAC ordering": 16
            },
            "Mn": {
                "Atomic mass": 54.938045,
                "Atomic no": 25,
                "Atomic orbitals": {
                    "1s": -233.696912,
                    "2p": -23.066297,
                    "2s": -26.866646,
                    "3d": -0.26654,
                    "3p": -1.99145,
                    "3s": -3.076637,
                    "4s": -0.191136
                },
                "Atomic radius": 1.4,
                "Atomic radius calculated": 1.61,
                "Boiling point": "2334 K",
                "Brinell hardness": "196 MN m<sup>-2</sup>",
                "Bulk modulus": "120 GPa",
                "Coefficient of linear thermal expansion": "21.7 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    2,
                    4,
                    7
                ],
                "Critical temperature": "no data K",
                "Density of solid": "7470 kg m<sup>-3</sup>",
                "Electrical resistivity": "144 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ar].3d<sup>5</sup>.4s<sup>2</sup>",
                "ICSD oxidation states": [
                    2,
                    3,
                    4,
                    7
                ],
                "Ionic radii": {
                    "2": 0.97,
                    "3": 0.785,
                    "4": 0.67,
                    "5": 0.47,
                    "6": 0.395,
                    "7": 0.6
                },
                "Ionic radii hs": {
                    "2": 0.97,
                    "3": 0.785
                },
                "Ionic radii ls": {
                    "2": 0.81,
                    "3": 0.72,
                    "4": 0.67,
                    "5": 0.47,
                    "6": 0.395,
                    "7": 0.6
                },
                "Liquid range": "815 K",
                "Melting point": "1519 K",
                "Mendeleev no": 60,
                "Mineral hardness": "6.0",
                "Molar volume": "7.35 cm<sup>3</sup>",
                "Name": "Manganese",
                "Oxidation states": [
                    -3,
                    -2,
                    -1,
                    1,
                    2,
                    3,
                    4,
                    5,
                    6,
                    7
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "2": {
                        "IV": {
                            "High Spin": {
                                "crystal_radius": 0.8,
                                "ionic_radius": 0.66
                            }
                        },
                        "V": {
                            "High Spin": {
                                "crystal_radius": 0.89,
                                "ionic_radius": 0.75
                            }
                        },
                        "VI": {
                            "Low Spin": {
                                "crystal_radius": 0.81,
                                "ionic_radius": 0.67
                            },
                            "High Spin": {
                                "crystal_radius": 0.97,
                                "ionic_radius": 0.83
                            }
                        },
                        "VII": {
                            "High Spin": {
                                "crystal_radius": 1.04,
                                "ionic_radius": 0.9
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.1,
                                "ionic_radius": 0.96
                            }
                        }
                    },
                    "3": {
                        "V": {
                            "": {
                                "crystal_radius": 0.72,
                                "ionic_radius": 0.58
                            }
                        },
                        "VI": {
                            "Low Spin": {
                                "crystal_radius": 0.72,
                                "ionic_radius": 0.58
                            },
                            "High Spin": {
                                "crystal_radius": 0.785,
                                "ionic_radius": 0.645
                            }
                        }
                    },
                    "4": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.53,
                                "ionic_radius": 0.39
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.67,
                                "ionic_radius": 0.53
                            }
                        }
                    },
                    "5": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.47,
                                "ionic_radius": 0.33
                            }
                        }
                    },
                    "6": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.395,
                                "ionic_radius": 0.255
                            }
                        }
                    },
                    "7": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.39,
                                "ionic_radius": 0.25
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.6,
                                "ionic_radius": 0.46
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "7.8 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "5150 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.55,
                "Youngs modulus": "198 GPa",
                "NMR Quadrupole Moment": {
                    "Mn-55": 330.1
                },
                "Metallic radius": 1.292,
                "iupac_ordering": 61,
                "IUPAC ordering": 61
            },
            "Mo": {
                "Atomic mass": 95.94,
                "Atomic no": 42,
                "Atomic orbitals": {
                    "1s": -709.232119,
                    "2p": -90.791541,
                    "2s": -98.503638,
                    "3d": -8.257721,
                    "3p": -13.71481,
                    "3s": -16.681545,
                    "4d": -0.153347,
                    "4p": -1.39005,
                    "4s": -2.234824,
                    "5s": -0.14788
                },
                "Atomic radius": 1.45,
                "Atomic radius calculated": 1.9,
                "Boiling point": "4912 K",
                "Brinell hardness": "1500 MN m<sup>-2</sup>",
                "Bulk modulus": "230 GPa",
                "Coefficient of linear thermal expansion": "4.8 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    4,
                    6
                ],
                "Critical temperature": "no data K",
                "Density of solid": "10280 kg m<sup>-3</sup>",
                "Electrical resistivity": "5.5 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Kr].4d<sup>5</sup>.5s<sup>1</sup>",
                "ICSD oxidation states": [
                    2,
                    3,
                    4,
                    5,
                    6
                ],
                "Ionic radii": {
                    "3": 0.83,
                    "4": 0.79,
                    "5": 0.75,
                    "6": 0.73
                },
                "Liquid range": "2016 K",
                "Melting point": "2896 K",
                "Mendeleev no": 56,
                "Mineral hardness": "5.5",
                "Molar volume": "9.38 cm<sup>3</sup>",
                "Name": "Molybdenum",
                "Oxidation states": [
                    -2,
                    -1,
                    1,
                    2,
                    3,
                    4,
                    5,
                    6
                ],
                "Poissons ratio": "0.31",
                "Reflectivity": "58 %",
                "Refractive index": "no data",
                "Rigidity modulus": "20 GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.83,
                                "ionic_radius": 0.69
                            }
                        }
                    },
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.79,
                                "ionic_radius": 0.65
                            }
                        }
                    },
                    "5": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.6,
                                "ionic_radius": 0.46
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.75,
                                "ionic_radius": 0.61
                            }
                        }
                    },
                    "6": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.55,
                                "ionic_radius": 0.41
                            }
                        },
                        "V": {
                            "": {
                                "crystal_radius": 0.64,
                                "ionic_radius": 0.5
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.73,
                                "ionic_radius": 0.59
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 0.87,
                                "ionic_radius": 0.73
                            }
                        }
                    }
                },
                "Superconduction temperature": "0.915 K",
                "Thermal conductivity": "139 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "6190 m s<sup>-1</sup>",
                "Vickers hardness": "1530 MN m<sup>-2</sup>",
                "X": 2.16,
                "Youngs modulus": "329 GPa",
                "Metallic radius": 1.402,
                "iupac_ordering": 57,
                "IUPAC ordering": 57
            },
            "N": {
                "Atomic mass": 14.0067,
                "Atomic no": 7,
                "Atomic orbitals": {
                    "1s": -14.011501,
                    "2p": -0.266297,
                    "2s": -0.676151
                },
                "Atomic radius": 0.65,
                "Atomic radius calculated": 0.56,
                "Boiling point": "77.36 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    -3,
                    3,
                    5
                ],
                "Critical temperature": "126.2 K",
                "Density of solid": "no data kg m<sup>-3</sup>",
                "Electrical resistivity": "no data 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[He].2s<sup>2</sup>.2p<sup>3</sup>",
                "ICSD oxidation states": [
                    1,
                    3,
                    5,
                    -1,
                    -3,
                    -2
                ],
                "Ionic radii": {
                    "-3": 1.32,
                    "3": 0.3,
                    "5": 0.27
                },
                "Liquid range": "14.31 K",
                "Melting point": "63.05 K",
                "Mendeleev no": 100,
                "Mineral hardness": "no data",
                "Molar volume": "13.54 cm<sup>3</sup>",
                "Name": "Nitrogen",
                "Oxidation states": [
                    -3,
                    -2,
                    -1,
                    1,
                    2,
                    3,
                    4,
                    5
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "1.000298 (gas; liquid 1.197)(no units)",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "-3": {
                        "IV": {
                            "": {
                                "crystal_radius": 1.32,
                                "ionic_radius": 1.46
                            }
                        }
                    },
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.3,
                                "ionic_radius": 0.16
                            }
                        }
                    },
                    "5": {
                        "III": {
                            "": {
                                "crystal_radius": 0.044,
                                "ionic_radius": -0.104
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.27,
                                "ionic_radius": 0.13
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "0.02583 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.55,
                "Velocity of sound": "333.6 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 3.04,
                "Youngs modulus": "no data GPa",
                "NMR Quadrupole Moment": {
                    "N-14": 20.44
                },
                "Metallic radius": "no data",
                "iupac_ordering": 91,
                "IUPAC ordering": 91
            },
            "Na": {
                "Atomic mass": 22.98976928,
                "Atomic no": 11,
                "Atomic orbitals": {
                    "1s": -37.719975,
                    "2p": -1.060636,
                    "2s": -2.063401,
                    "3s": -0.103415
                },
                "Atomic radius": 1.8,
                "Atomic radius calculated": 1.9,
                "Boiling point": "1156 K",
                "Brinell hardness": "0.69 MN m<sup>-2</sup>",
                "Bulk modulus": "6.3 GPa",
                "Coefficient of linear thermal expansion": "71 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    1
                ],
                "Critical temperature": "2573 K",
                "Density of solid": "968 kg m<sup>-3</sup>",
                "Electrical resistivity": "4.9 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ne].3s<sup>1</sup>",
                "ICSD oxidation states": [
                    1
                ],
                "Ionic radii": {
                    "1": 1.16
                },
                "Liquid range": "785.13 K",
                "Melting point": "370.87 K",
                "Mendeleev no": 11,
                "Mineral hardness": "0.5",
                "Molar volume": "23.78 cm<sup>3</sup>",
                "Name": "Sodium",
                "Oxidation states": [
                    -1,
                    1
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "3.3 GPa",
                "Shannon radii": {
                    "1": {
                        "IV": {
                            "": {
                                "crystal_radius": 1.13,
                                "ionic_radius": 0.99
                            }
                        },
                        "V": {
                            "": {
                                "crystal_radius": 1.14,
                                "ionic_radius": 1.0
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 1.16,
                                "ionic_radius": 1.02
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.26,
                                "ionic_radius": 1.12
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.32,
                                "ionic_radius": 1.18
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.38,
                                "ionic_radius": 1.24
                            }
                        },
                        "XII": {
                            "": {
                                "crystal_radius": 1.53,
                                "ionic_radius": 1.39
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "140 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 2.27,
                "Velocity of sound": "3200 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 0.93,
                "Youngs modulus": "10 GPa",
                "NMR Quadrupole Moment": {
                    "Na-23": 104.1
                },
                "Metallic radius": 1.86,
                "iupac_ordering": 10,
                "IUPAC ordering": 10
            },
            "Nb": {
                "Atomic mass": 92.90638,
                "Atomic no": 41,
                "Atomic orbitals": {
                    "1s": -673.76253,
                    "2p": -85.272175,
                    "2s": -92.74086,
                    "3d": -7.339839,
                    "3p": -12.552855,
                    "3s": -15.393727,
                    "4d": -0.125252,
                    "4p": -1.250049,
                    "4s": -2.036693,
                    "5s": -0.144272
                },
                "Atomic radius": 1.45,
                "Atomic radius calculated": 1.98,
                "Boiling point": "5017 K",
                "Brinell hardness": "736 MN m<sup>-2</sup>",
                "Bulk modulus": "170 GPa",
                "Coefficient of linear thermal expansion": "7.3 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    5
                ],
                "Critical temperature": "no data K",
                "Density of solid": "8570 kg m<sup>-3</sup>",
                "Electrical resistivity": "15.2 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Kr].4d<sup>4</sup>.5s<sup>1</sup>",
                "ICSD oxidation states": [
                    2,
                    3,
                    4,
                    5
                ],
                "Ionic radii": {
                    "3": 0.86,
                    "4": 0.82,
                    "5": 0.78
                },
                "Liquid range": "2267 K",
                "Melting point": "2750 K",
                "Mendeleev no": 53,
                "Mineral hardness": "6.0",
                "Molar volume": "10.83 cm<sup>3</sup>",
                "Name": "Niobium",
                "Oxidation states": [
                    -1,
                    2,
                    3,
                    4,
                    5
                ],
                "Poissons ratio": "0.40",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "38 GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.86,
                                "ionic_radius": 0.72
                            }
                        }
                    },
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.82,
                                "ionic_radius": 0.68
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 0.93,
                                "ionic_radius": 0.79
                            }
                        }
                    },
                    "5": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.62,
                                "ionic_radius": 0.48
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.78,
                                "ionic_radius": 0.64
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 0.83,
                                "ionic_radius": 0.69
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 0.88,
                                "ionic_radius": 0.74
                            }
                        }
                    }
                },
                "Superconduction temperature": "9.25 K",
                "Thermal conductivity": "54 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "3480 m s<sup>-1</sup>",
                "Vickers hardness": "1320 MN m<sup>-2</sup>",
                "X": 1.6,
                "Youngs modulus": "105 GPa",
                "Metallic radius": 1.473,
                "iupac_ordering": 54,
                "IUPAC ordering": 54
            },
            "Nd": {
                "Atomic mass": 144.242,
                "Atomic no": 60,
                "Atomic orbitals": {
                    "1s": -1509.698955,
                    "2p": -224.351816,
                    "2s": -236.613572,
                    "3d": -35.754515,
                    "3p": -45.791219,
                    "3s": -51.161263,
                    "4d": -4.377027,
                    "4f": -0.179508,
                    "4p": -7.96782,
                    "4s": -10.000891,
                    "5p": -0.798503,
                    "5s": -1.334934,
                    "6s": -0.125796
                },
                "Atomic radius": 1.85,
                "Atomic radius calculated": 2.06,
                "Boiling point": "3373 K",
                "Brinell hardness": "265 MN m<sup>-2</sup>",
                "Bulk modulus": "32 GPa",
                "Coefficient of linear thermal expansion": "9.6 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "6800 kg m<sup>-3</sup>",
                "Electrical resistivity": "64.3 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>4</sup>.6s<sup>2</sup>",
                "ICSD oxidation states": [
                    2,
                    3
                ],
                "Ionic radii": {
                    "2": 1.43,
                    "3": 1.123
                },
                "Liquid range": "2076 K",
                "Melting point": "1297 K",
                "Mendeleev no": 30,
                "Mineral hardness": "no data",
                "Molar volume": "20.59 cm<sup>3</sup>",
                "Name": "Neodymium",
                "Oxidation states": [
                    2,
                    3
                ],
                "Poissons ratio": "0.28",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "16 GPa",
                "Shannon radii": {
                    "2": {
                        "VIII": {
                            "": {
                                "crystal_radius": 1.43,
                                "ionic_radius": 1.29
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.49,
                                "ionic_radius": 1.35
                            }
                        }
                    },
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.123,
                                "ionic_radius": 0.983
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.249,
                                "ionic_radius": 1.109
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.303,
                                "ionic_radius": 1.163
                            }
                        },
                        "XII": {
                            "": {
                                "crystal_radius": 1.41,
                                "ionic_radius": 1.27
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "17 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "2330 m s<sup>-1</sup>",
                "Vickers hardness": "343 MN m<sup>-2</sup>",
                "X": 1.14,
                "Youngs modulus": "41 GPa",
                "Metallic radius": 1.821,
                "iupac_ordering": 44,
                "IUPAC ordering": 44
            },
            "Ne": {
                "Atomic mass": 20.1797,
                "Atomic no": 10,
                "Atomic orbitals": {
                    "1s": -30.305855,
                    "2p": -0.498034,
                    "2s": -1.322809
                },
                "Atomic radius": "no data",
                "Atomic radius calculated": 0.38,
                "Boiling point": "27.07 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Critical temperature": "44.4 K",
                "Density of solid": "no data kg m<sup>-3</sup>",
                "Electrical resistivity": "no data 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[He].2s<sup>2</sup>.2p<sup>6</sup>",
                "Liquid range": "2.51 K",
                "Max oxidation state": 0.0,
                "Melting point": "24.56 K",
                "Mendeleev no": 2,
                "Min oxidation state": 0.0,
                "Mineral hardness": "no data",
                "Molar volume": "13.23 cm<sup>3</sup>",
                "Name": "Neon",
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "1.000067",
                "Rigidity modulus": "no data GPa",
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "0.0491 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.54,
                "Velocity of sound": "936 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "Youngs modulus": "no data GPa",
                "NMR Quadrupole Moment": {
                    "Ne-21": 101.55
                },
                "Metallic radius": "no data",
                "iupac_ordering": 4,
                "IUPAC ordering": 4
            },
            "Ni": {
                "Atomic mass": 58.6934,
                "Atomic no": 28,
                "Atomic orbitals": {
                    "1s": -297.870824,
                    "2p": -30.868027,
                    "2s": -35.312112,
                    "3d": -0.348699,
                    "3p": -2.594158,
                    "3s": -3.950717,
                    "4s": -0.210764
                },
                "Atomic radius": 1.35,
                "Atomic radius calculated": 1.49,
                "Boiling point": "3186 K",
                "Brinell hardness": "700 MN m<sup>-2</sup>",
                "Bulk modulus": "180 GPa",
                "Coefficient of linear thermal expansion": "13.4 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    2
                ],
                "Critical temperature": "no data K",
                "Density of solid": "8908 kg m<sup>-3</sup>",
                "Electrical resistivity": "7.2 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ar].3d<sup>8</sup>.4s<sup>2</sup>",
                "ICSD oxidation states": [
                    1,
                    2,
                    3,
                    4
                ],
                "Ionic radii": {
                    "3": 0.74
                },
                "Ionic radii hs": {
                    "3": 0.74
                },
                "Ionic radii ls": {
                    "2": 0.83,
                    "3": 0.7,
                    "4": 0.62
                },
                "Liquid range": "1458 K",
                "Melting point": "1728 K",
                "Mendeleev no": 67,
                "Mineral hardness": "4.0",
                "Molar volume": "6.59 cm<sup>3</sup>",
                "Name": "Nickel",
                "Oxidation states": [
                    -1,
                    1,
                    2,
                    3,
                    4
                ],
                "Poissons ratio": "0.31",
                "Reflectivity": "72 %",
                "Refractive index": "no data",
                "Rigidity modulus": "76 GPa",
                "Shannon radii": {
                    "2": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.69,
                                "ionic_radius": 0.55
                            }
                        },
                        "IVSQ": {
                            "": {
                                "crystal_radius": 0.63,
                                "ionic_radius": 0.49
                            }
                        },
                        "V": {
                            "": {
                                "crystal_radius": 0.77,
                                "ionic_radius": 0.63
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.83,
                                "ionic_radius": 0.69
                            }
                        }
                    },
                    "3": {
                        "VI": {
                            "Low Spin": {
                                "crystal_radius": 0.7,
                                "ionic_radius": 0.56
                            },
                            "High Spin": {
                                "crystal_radius": 0.74,
                                "ionic_radius": 0.6
                            }
                        }
                    },
                    "4": {
                        "VI": {
                            "Low Spin": {
                                "crystal_radius": 0.62,
                                "ionic_radius": 0.48
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "91 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.63,
                "Velocity of sound": "4970 m s<sup>-1</sup>",
                "Vickers hardness": "638 MN m<sup>-2</sup>",
                "X": 1.91,
                "Youngs modulus": "200 GPa",
                "NMR Quadrupole Moment": {
                    "Ni-61": 162.15
                },
                "Metallic radius": 1.246,
                "iupac_ordering": 70,
                "IUPAC ordering": 70
            },
            "No": {
                "Atomic mass": 259.0,
                "Atomic no": 102,
                "Atomic orbitals": "no data",
                "Atomic radius": "no data",
                "Atomic radius calculated": "no data",
                "Boiling point": "no data K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "no data kg m<sup>-3</sup>",
                "Electrical resistivity": "no data 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Rn].5f<sup>14</sup>.7s<sup>2</sup>",
                "Liquid range": "no data K",
                "Melting point": "about 1100 K",
                "Mendeleev no": 35,
                "Mineral hardness": "no data",
                "Molar volume": "no data cm<sup>3</sup>",
                "Name": "Nobelium",
                "Oxidation states": [
                    2,
                    3
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "2": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.24,
                                "ionic_radius": 1.1
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "no data W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.3,
                "Youngs modulus": "no data GPa",
                "Metallic radius": "no data",
                "iupac_ordering": 19,
                "IUPAC ordering": 19
            },
            "Np": {
                "Atomic mass": 237.0,
                "Atomic no": 93,
                "Atomic orbitals": "no data",
                "Atomic radius": 1.75,
                "Atomic radius calculated": "no data",
                "Boiling point": "4273 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    5
                ],
                "Critical temperature": "no data K",
                "Density of solid": "20450 kg m<sup>-3</sup>",
                "Electrical resistivity": "120 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Rn].5f<sup>4</sup>.6d<sup>1</sup>.7s<sup>2</sup>",
                "Ionic radii": {
                    "2": 1.24,
                    "3": 1.15,
                    "4": 1.01,
                    "5": 0.89,
                    "6": 0.86,
                    "7": 0.85
                },
                "Liquid range": "3363 K",
                "Melting point": "910 K",
                "Mendeleev no": 44,
                "Mineral hardness": "no data",
                "Molar volume": "11.59 cm<sup>3</sup>",
                "Name": "Neptunium",
                "Oxidation states": [
                    3,
                    4,
                    5,
                    6,
                    7
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "2": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.24,
                                "ionic_radius": 1.1
                            }
                        }
                    },
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.15,
                                "ionic_radius": 1.01
                            }
                        }
                    },
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.01,
                                "ionic_radius": 0.87
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.12,
                                "ionic_radius": 0.98
                            }
                        }
                    },
                    "5": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.89,
                                "ionic_radius": 0.75
                            }
                        }
                    },
                    "6": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.86,
                                "ionic_radius": 0.72
                            }
                        }
                    },
                    "7": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.85,
                                "ionic_radius": 0.71
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "6 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.36,
                "Youngs modulus": "no data GPa",
                "Metallic radius": 1.503,
                "iupac_ordering": 28,
                "IUPAC ordering": 28
            },
            "O": {
                "Atomic mass": 15.9994,
                "Atomic no": 8,
                "Atomic orbitals": {
                    "1s": -18.758245,
                    "2p": -0.338381,
                    "2s": -0.871362
                },
                "Atomic radius": 0.6,
                "Atomic radius calculated": 0.48,
                "Boiling point": "90.2 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    -2
                ],
                "Critical temperature": "154.6 K",
                "Density of solid": "no data kg m<sup>-3</sup>",
                "Electrical resistivity": "no data 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[He].2s<sup>2</sup>.2p<sup>4</sup>",
                "ICSD oxidation states": [
                    -2
                ],
                "Ionic radii": {
                    "-2": 1.26
                },
                "Liquid range": "35.4 K",
                "Melting point": "54.8 K",
                "Mendeleev no": 101,
                "Mineral hardness": "no data",
                "Molar volume": "17.36 cm<sup>3</sup>",
                "Name": "Oxygen",
                "Oxidation states": [
                    -2,
                    -1,
                    1,
                    2
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "1.000271 (gas; liquid 1.221)(no units)",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "-2": {
                        "II": {
                            "": {
                                "crystal_radius": 1.21,
                                "ionic_radius": 1.35
                            }
                        },
                        "III": {
                            "": {
                                "crystal_radius": 1.22,
                                "ionic_radius": 1.36
                            }
                        },
                        "IV": {
                            "": {
                                "crystal_radius": 1.24,
                                "ionic_radius": 1.38
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 1.26,
                                "ionic_radius": 1.4
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.28,
                                "ionic_radius": 1.42
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "0.02658 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.52,
                "Velocity of sound": "317.5 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 3.44,
                "Youngs modulus": "no data GPa",
                "NMR Quadrupole Moment": {
                    "O-17": -25.58
                },
                "Metallic radius": "no data",
                "iupac_ordering": 97,
                "IUPAC ordering": 97
            },
            "Os": {
                "Atomic mass": 190.23,
                "Atomic no": 76,
                "Atomic orbitals": {
                    "1s": -2475.238617,
                    "2p": -393.15408,
                    "2s": -409.522396,
                    "3d": -72.497183,
                    "3p": -86.837047,
                    "3s": -94.501324,
                    "4d": -10.176082,
                    "4f": -2.321175,
                    "4p": -16.119671,
                    "4s": -19.362527,
                    "5d": -0.296791,
                    "5p": -1.757404,
                    "5s": -2.738293,
                    "6s": -0.191489
                },
                "Atomic radius": 1.3,
                "Atomic radius calculated": 1.85,
                "Boiling point": "5285 K",
                "Brinell hardness": "3920 MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "5.1 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    4
                ],
                "Critical temperature": "no data K",
                "Density of solid": "22610 kg m<sup>-3</sup>",
                "Electrical resistivity": "8.1 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>14</sup>.5d<sup>6</sup>.6s<sup>2</sup>",
                "Ionic radii": {
                    "4": 0.77,
                    "5": 0.715,
                    "6": 0.685,
                    "7": 0.665,
                    "8": 0.53
                },
                "Liquid range": "1979 K",
                "Melting point": "3306 K",
                "Mendeleev no": 63,
                "Mineral hardness": "7.0",
                "Molar volume": "8.42 cm<sup>3</sup>",
                "Name": "Osmium",
                "Oxidation states": [
                    -2,
                    -1,
                    1,
                    2,
                    3,
                    4,
                    5,
                    6,
                    7,
                    8
                ],
                "Poissons ratio": "0.25",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "222 GPa",
                "Shannon radii": {
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.77,
                                "ionic_radius": 0.63
                            }
                        }
                    },
                    "5": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.715,
                                "ionic_radius": 0.575
                            }
                        }
                    },
                    "6": {
                        "V": {
                            "": {
                                "crystal_radius": 0.63,
                                "ionic_radius": 0.49
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.685,
                                "ionic_radius": 0.545
                            }
                        }
                    },
                    "7": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.665,
                                "ionic_radius": 0.525
                            }
                        }
                    },
                    "8": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.53,
                                "ionic_radius": 0.39
                            }
                        }
                    }
                },
                "Superconduction temperature": "0.66 K",
                "Thermal conductivity": "88 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "4940 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 2.2,
                "Youngs modulus": "no data GPa",
                "Metallic radius": 1.352,
                "iupac_ordering": 62,
                "IUPAC ordering": 62
            },
            "P": {
                "Atomic mass": 30.973762,
                "Atomic no": 15,
                "Atomic orbitals": {
                    "1s": -76.061897,
                    "2p": -4.576617,
                    "2s": -6.329346,
                    "3p": -0.20608,
                    "3s": -0.512364
                },
                "Atomic radius": 1.0,
                "Atomic radius calculated": 0.98,
                "Boiling point": "550 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "11 GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    -3,
                    3,
                    5
                ],
                "Critical temperature": "994 K",
                "Density of solid": "1823 kg m<sup>-3</sup>",
                "Electrical resistivity": "about 10 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ne].3s<sup>2</sup>.3p<sup>3</sup>",
                "ICSD oxidation states": [
                    3,
                    4,
                    5,
                    -2,
                    -3,
                    -1
                ],
                "Ionic radii": {
                    "3": 0.58,
                    "5": 0.52
                },
                "Liquid range": "232.7 K",
                "Melting point": "(white P) 317.3 K",
                "Mendeleev no": 90,
                "Mineral hardness": "no data",
                "Molar volume": "17.02 cm<sup>3</sup>",
                "Name": "Phosphorus",
                "Oxidation states": [
                    -3,
                    -2,
                    -1,
                    1,
                    2,
                    3,
                    4,
                    5
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "1.001212",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.58,
                                "ionic_radius": 0.44
                            }
                        }
                    },
                    "5": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.31,
                                "ionic_radius": 0.17
                            }
                        },
                        "V": {
                            "": {
                                "crystal_radius": 0.43,
                                "ionic_radius": 0.29
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.52,
                                "ionic_radius": 0.38
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "0.236 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.8,
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 2.19,
                "Youngs modulus": "no data GPa",
                "Metallic radius": "no data",
                "iupac_ordering": 90,
                "IUPAC ordering": 90
            },
            "Pa": {
                "Atomic mass": 231.03588,
                "Atomic no": 91,
                "Atomic orbitals": {
                    "1s": -3606.333629,
                    "2p": -603.470278,
                    "2s": -623.870431,
                    "3d": -127.781168,
                    "3p": -146.485678,
                    "3s": -156.466742,
                    "4d": -25.933121,
                    "4f": -14.105747,
                    "4p": -34.48293,
                    "4s": -39.064507,
                    "5d": -3.659928,
                    "5f": -0.316813,
                    "5p": -6.709821,
                    "5s": -8.463463,
                    "6d": -0.142481,
                    "6p": -0.799756,
                    "6s": -1.287232,
                    "7s": -0.129653
                },
                "Atomic radius": 1.8,
                "Atomic radius calculated": "no data",
                "Boiling point": "no data K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    5
                ],
                "Critical temperature": "no data K",
                "Density of solid": "15370 kg m<sup>-3</sup>",
                "Electrical resistivity": "18 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Rn].5f<sup>2</sup>.6d<sup>1</sup>.7s<sup>2</sup>",
                "Ionic radii": {
                    "3": 1.16,
                    "4": 1.04,
                    "5": 0.92
                },
                "Liquid range": "no data K",
                "Melting point": "1841 K",
                "Mendeleev no": 46,
                "Mineral hardness": "no data",
                "Molar volume": "15.18 cm<sup>3</sup>",
                "Name": "Protactinium",
                "Oxidation states": [
                    3,
                    4,
                    5
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.18,
                                "ionic_radius": 1.04
                            }
                        }
                    },
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.04,
                                "ionic_radius": 0.9
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.15,
                                "ionic_radius": 1.01
                            }
                        }
                    },
                    "5": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.92,
                                "ionic_radius": 0.78
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.05,
                                "ionic_radius": 0.91
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.09,
                                "ionic_radius": 0.95
                            }
                        }
                    }
                },
                "Superconduction temperature": "1.4 K",
                "Thermal conductivity": "47 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.5,
                "Youngs modulus": "no data GPa",
                "Metallic radius": 1.642,
                "iupac_ordering": 30,
                "IUPAC ordering": 30
            },
            "Pb": {
                "Atomic mass": 207.2,
                "Atomic no": 82,
                "Atomic orbitals": {
                    "1s": -2901.078061,
                    "2p": -470.877785,
                    "2s": -488.843335,
                    "3d": -91.889924,
                    "3p": -107.950391,
                    "3s": -116.526852,
                    "4d": -15.030026,
                    "4f": -5.592532,
                    "4p": -21.990564,
                    "4s": -25.75333,
                    "5d": -0.902393,
                    "5p": -2.941657,
                    "5s": -4.206797,
                    "6p": -0.141831,
                    "6s": -0.357187
                },
                "Atomic radius": 1.8,
                "Atomic radius calculated": 1.54,
                "Boiling point": "2022 K",
                "Brinell hardness": "38.3 MN m<sup>-2</sup>",
                "Bulk modulus": "46 GPa",
                "Coefficient of linear thermal expansion": "28.9 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    2,
                    4
                ],
                "Critical temperature": "no data K",
                "Density of solid": "11340 kg m<sup>-3</sup>",
                "Electrical resistivity": "21 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>14</sup>.5d<sup>10</sup>.6s<sup>2</sup>.6p<sup>2</sup>",
                "ICSD oxidation states": [
                    2,
                    4
                ],
                "Ionic radii": {
                    "2": 1.33,
                    "4": 0.915
                },
                "Liquid range": "1421.39 K",
                "Melting point": "600.61 K",
                "Mendeleev no": 82,
                "Mineral hardness": "1.5",
                "Molar volume": "18.26 cm<sup>3</sup>",
                "Name": "Lead",
                "Oxidation states": [
                    -4,
                    2,
                    4
                ],
                "Poissons ratio": "0.44",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "5.6 GPa",
                "Shannon radii": {
                    "2": {
                        "IVPY": {
                            "": {
                                "crystal_radius": 1.12,
                                "ionic_radius": 0.98
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 1.33,
                                "ionic_radius": 1.19
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.37,
                                "ionic_radius": 1.23
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.43,
                                "ionic_radius": 1.29
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.49,
                                "ionic_radius": 1.35
                            }
                        },
                        "X": {
                            "": {
                                "crystal_radius": 1.54,
                                "ionic_radius": 1.4
                            }
                        },
                        "XI": {
                            "": {
                                "crystal_radius": 1.59,
                                "ionic_radius": 1.45
                            }
                        },
                        "XII": {
                            "": {
                                "crystal_radius": 1.63,
                                "ionic_radius": 1.49
                            }
                        }
                    },
                    "4": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.79,
                                "ionic_radius": 0.65
                            }
                        },
                        "V": {
                            "": {
                                "crystal_radius": 0.87,
                                "ionic_radius": 0.73
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.915,
                                "ionic_radius": 0.775
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.08,
                                "ionic_radius": 0.94
                            }
                        }
                    }
                },
                "Superconduction temperature": "7.2 K",
                "Thermal conductivity": "35 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 2.02,
                "Velocity of sound": "1260 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 2.33,
                "Youngs modulus": "16 GPa",
                "Metallic radius": 1.75,
                "iupac_ordering": 82,
                "IUPAC ordering": 82
            },
            "Pd": {
                "Atomic mass": 106.42,
                "Atomic no": 46,
                "Atomic orbitals": {
                    "1s": -860.134909,
                    "2p": -114.408286,
                    "2s": -123.105078,
                    "3d": -12.132197,
                    "3p": -18.580798,
                    "3s": -22.060898,
                    "4d": -0.160771,
                    "4p": -1.815215,
                    "4s": -2.889173
                },
                "Atomic radius": 1.4,
                "Atomic radius calculated": 1.69,
                "Boiling point": "3236 K",
                "Brinell hardness": "37.3 MN m<sup>-2</sup>",
                "Bulk modulus": "180 GPa",
                "Coefficient of linear thermal expansion": "11.8 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    2,
                    4
                ],
                "Critical temperature": "no data K",
                "Density of solid": "12023 kg m<sup>-3</sup>",
                "Electrical resistivity": "10.8 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Kr].4d<sup>10</sup>",
                "ICSD oxidation states": [
                    2,
                    4
                ],
                "Ionic radii": {
                    "1": 0.73,
                    "2": 1.0,
                    "3": 0.9,
                    "4": 0.755
                },
                "Liquid range": "1407.95 K",
                "Melting point": "1828.05 K",
                "Mendeleev no": 69,
                "Mineral hardness": "4.75",
                "Molar volume": "8.56 cm<sup>3</sup>",
                "Name": "Palladium",
                "Oxidation states": [
                    2,
                    4
                ],
                "Poissons ratio": "0.39",
                "Reflectivity": "72 %",
                "Refractive index": "no data",
                "Rigidity modulus": "44 GPa",
                "Shannon radii": {
                    "1": {
                        "II": {
                            "": {
                                "crystal_radius": 0.73,
                                "ionic_radius": 0.59
                            }
                        }
                    },
                    "2": {
                        "IVSQ": {
                            "": {
                                "crystal_radius": 0.78,
                                "ionic_radius": 0.64
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 1.0,
                                "ionic_radius": 0.86
                            }
                        }
                    },
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.9,
                                "ionic_radius": 0.76
                            }
                        }
                    },
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.755,
                                "ionic_radius": 0.615
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "72 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.63,
                "Velocity of sound": "3070 m s<sup>-1</sup>",
                "Vickers hardness": "461 MN m<sup>-2</sup>",
                "X": 2.2,
                "Youngs modulus": "121 GPa",
                "Metallic radius": 1.376,
                "iupac_ordering": 69,
                "IUPAC ordering": 69
            },
            "Pm": {
                "Atomic mass": 145.0,
                "Atomic no": 61,
                "Atomic orbitals": {
                    "1s": -1562.980284,
                    "2p": -233.455114,
                    "2s": -245.970548,
                    "3d": -37.625433,
                    "3p": -47.921132,
                    "3s": -53.429311,
                    "4d": -4.596822,
                    "4f": -0.200159,
                    "4p": -8.320495,
                    "4s": -10.422756,
                    "5p": -0.817702,
                    "5s": -1.372265,
                    "6s": -0.127053
                },
                "Atomic radius": 1.85,
                "Atomic radius calculated": 2.05,
                "Boiling point": "3273 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "33 GPa",
                "Coefficient of linear thermal expansion": "11 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "7264 kg m<sup>-3</sup>",
                "Electrical resistivity": "about 75 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>5</sup>.6s<sup>2</sup>",
                "Ionic radii": {
                    "3": 1.11
                },
                "Liquid range": "1900 K",
                "Melting point": "1373 K",
                "Mendeleev no": 29,
                "Mineral hardness": "no data",
                "Molar volume": "20.23 cm<sup>3</sup>",
                "Name": "Promethium",
                "Oxidation states": [
                    3
                ],
                "Poissons ratio": "0.28",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "18 GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.11,
                                "ionic_radius": 0.97
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.233,
                                "ionic_radius": 1.093
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.284,
                                "ionic_radius": 1.144
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "15 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.13,
                "Youngs modulus": "46 GPa",
                "Metallic radius": 1.811,
                "iupac_ordering": 43,
                "IUPAC ordering": 43
            },
            "Po": {
                "Atomic mass": 210.0,
                "Atomic no": 84,
                "Atomic orbitals": {
                    "1s": -3050.988417,
                    "2p": -498.77192,
                    "2s": -517.275843,
                    "3d": -99.256068,
                    "3p": -115.898384,
                    "3s": -124.783683,
                    "4d": -17.173307,
                    "4f": -7.206499,
                    "4p": -24.481337,
                    "4s": -28.42254,
                    "5d": -1.386458,
                    "5p": -3.655382,
                    "5s": -5.027447,
                    "6p": -0.217889,
                    "6s": -0.493528
                },
                "Atomic radius": 1.9,
                "Atomic radius calculated": 1.35,
                "Boiling point": "1235 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    -2,
                    2,
                    4
                ],
                "Critical temperature": "no data K",
                "Density of solid": "9196 kg m<sup>-3</sup>",
                "Electrical resistivity": "40 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>14</sup>.5d<sup>10</sup>.6s<sup>2</sup>.6p<sup>4</sup>",
                "Ionic radii": {
                    "4": 1.08,
                    "6": 0.81
                },
                "Liquid range": "708 K",
                "Melting point": "527 K",
                "Mendeleev no": 91,
                "Mineral hardness": "no data",
                "Molar volume": "22.97 cm<sup>3</sup>",
                "Name": "Polonium",
                "Oxidation states": [
                    -2,
                    2,
                    4,
                    6
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.08,
                                "ionic_radius": 0.94
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.22,
                                "ionic_radius": 1.08
                            }
                        }
                    },
                    "6": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.81,
                                "ionic_radius": 0.67
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "20 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.97,
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 2.0,
                "Youngs modulus": "no data GPa",
                "Metallic radius": 1.53,
                "iupac_ordering": 93,
                "IUPAC ordering": 93
            },
            "Pr": {
                "Atomic mass": 140.90765,
                "Atomic no": 59,
                "Atomic orbitals": {
                    "1s": -1457.338067,
                    "2p": -215.418313,
                    "2s": -227.426363,
                    "3d": -33.913996,
                    "3p": -43.692548,
                    "3s": -48.924994,
                    "4d": -4.154228,
                    "4f": -0.155138,
                    "4p": -7.613108,
                    "4s": -9.577447,
                    "5p": -0.778046,
                    "5s": -1.296106,
                    "6s": -0.124465
                },
                "Atomic radius": 1.85,
                "Atomic radius calculated": 2.47,
                "Boiling point": "3563 K",
                "Brinell hardness": "481 MN m<sup>-2</sup>",
                "Bulk modulus": "29 GPa",
                "Coefficient of linear thermal expansion": "6.7 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "6640 kg m<sup>-3</sup>",
                "Electrical resistivity": "70 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>3</sup>.6s<sup>2</sup>",
                "ICSD oxidation states": [
                    3,
                    4
                ],
                "Ionic radii": {
                    "3": 1.13,
                    "4": 0.99
                },
                "Liquid range": "2355 K",
                "Melting point": "1208 K",
                "Mendeleev no": 31,
                "Mineral hardness": "no data",
                "Molar volume": "20.80 cm<sup>3</sup>",
                "Name": "Praseodymium",
                "Oxidation states": [
                    2,
                    3,
                    4
                ],
                "Poissons ratio": "0.28",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "15 GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.13,
                                "ionic_radius": 0.99
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.266,
                                "ionic_radius": 1.126
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.319,
                                "ionic_radius": 1.179
                            }
                        }
                    },
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.99,
                                "ionic_radius": 0.85
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.1,
                                "ionic_radius": 0.96
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "13 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "2280 m s<sup>-1</sup>",
                "Vickers hardness": "400 MN m<sup>-2</sup>",
                "X": 1.13,
                "Youngs modulus": "37 GPa",
                "Metallic radius": 1.828,
                "iupac_ordering": 45,
                "IUPAC ordering": 45
            },
            "Pt": {
                "Atomic mass": 195.084,
                "Atomic no": 78,
                "Atomic orbitals": {
                    "1s": -2613.096532,
                    "2p": -417.96053,
                    "2s": -434.858003,
                    "3d": -78.400271,
                    "3p": -93.309108,
                    "3s": -101.274869,
                    "4d": -11.419476,
                    "4f": -3.038049,
                    "4p": -17.697297,
                    "4s": -21.110651,
                    "5d": -0.273634,
                    "5p": -1.884256,
                    "5s": -2.950526,
                    "6s": -0.161308
                },
                "Atomic radius": 1.35,
                "Atomic radius calculated": 1.77,
                "Boiling point": "4098 K",
                "Brinell hardness": "392 MN m<sup>-2</sup>",
                "Bulk modulus": "230 GPa",
                "Coefficient of linear thermal expansion": "8.8 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    2,
                    4
                ],
                "Critical temperature": "no data K",
                "Density of solid": "21090 kg m<sup>-3</sup>",
                "Electrical resistivity": "10.6 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>14</sup>.5d<sup>9</sup>.6s<sup>1</sup>",
                "Ionic radii": {
                    "2": 0.94,
                    "4": 0.765,
                    "5": 0.71
                },
                "Liquid range": "2056.6 K",
                "Melting point": "2041.4 K",
                "Mendeleev no": 68,
                "Mineral hardness": "3.5",
                "Molar volume": "9.09 cm<sup>3</sup>",
                "Name": "Platinum",
                "Oxidation states": [
                    -2,
                    2,
                    4,
                    5,
                    6
                ],
                "Poissons ratio": "0.38",
                "Reflectivity": "73 %",
                "Refractive index": "no data",
                "Rigidity modulus": "61 GPa",
                "Shannon radii": {
                    "2": {
                        "IVSQ": {
                            "": {
                                "crystal_radius": 0.74,
                                "ionic_radius": 0.6
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.94,
                                "ionic_radius": 0.8
                            }
                        }
                    },
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.765,
                                "ionic_radius": 0.625
                            }
                        }
                    },
                    "5": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.71,
                                "ionic_radius": 0.57
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "72 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.75,
                "Velocity of sound": "2680 m s<sup>-1</sup>",
                "Vickers hardness": "549 MN m<sup>-2</sup>",
                "X": 2.28,
                "Youngs modulus": "168 GPa",
                "Metallic radius": 1.387,
                "iupac_ordering": 68,
                "IUPAC ordering": 68
            },
            "Pu": {
                "Atomic mass": 244.0,
                "Atomic no": 94,
                "Atomic orbitals": "no data",
                "Atomic radius": 1.75,
                "Atomic radius calculated": "no data",
                "Boiling point": "3503 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    4
                ],
                "Critical temperature": "no data K",
                "Density of solid": "19816 kg m<sup>-3</sup>",
                "Electrical resistivity": "150 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Rn].5f<sup>6</sup>.7s<sup>2</sup>",
                "Ionic radii": {
                    "3": 1.14,
                    "4": 1.0,
                    "5": 0.88,
                    "6": 0.85
                },
                "Liquid range": "2590.5 K",
                "Melting point": "912.5 K",
                "Mendeleev no": 43,
                "Mineral hardness": "no data",
                "Molar volume": "12.29 cm<sup>3</sup>",
                "Name": "Plutonium",
                "Oxidation states": [
                    3,
                    4,
                    5,
                    6,
                    7
                ],
                "Poissons ratio": "0.21",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "43 GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.14,
                                "ionic_radius": 1.0
                            }
                        }
                    },
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.0,
                                "ionic_radius": 0.86
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.1,
                                "ionic_radius": 0.96
                            }
                        }
                    },
                    "5": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.88,
                                "ionic_radius": 0.74
                            }
                        }
                    },
                    "6": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.85,
                                "ionic_radius": 0.71
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "6 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "2260 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.28,
                "Youngs modulus": "96 GPa",
                "Metallic radius": 1.523,
                "iupac_ordering": 27,
                "IUPAC ordering": 27
            },
            "Ra": {
                "Atomic mass": 226.0,
                "Atomic no": 88,
                "Atomic orbitals": {
                    "1s": -3362.736563,
                    "2p": -557.513214,
                    "2s": -577.101208,
                    "3d": -115.306476,
                    "3p": -133.12325,
                    "3s": -142.632426,
                    "4d": -22.208125,
                    "4f": -11.181066,
                    "4p": -30.221208,
                    "4s": -34.525628,
                    "5d": -2.819853,
                    "5p": -5.547203,
                    "5s": -7.139137,
                    "6p": -0.634674,
                    "6s": -1.05135,
                    "7s": -0.113732
                },
                "Atomic radius": 2.15,
                "Atomic radius calculated": "no data",
                "Boiling point": "2010 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    2
                ],
                "Critical temperature": "no data K",
                "Density of solid": "5000 kg m<sup>-3</sup>",
                "Electrical resistivity": "100 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Rn].7s<sup>2</sup>",
                "Ionic radii": {
                    "2": 1.62
                },
                "Liquid range": "1037 K",
                "Melting point": "973 K",
                "Mendeleev no": 13,
                "Mineral hardness": "no data",
                "Molar volume": "41.09 cm<sup>3</sup>",
                "Name": "Radium",
                "Oxidation states": [
                    2
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "2": {
                        "VIII": {
                            "": {
                                "crystal_radius": 1.62,
                                "ionic_radius": 1.48
                            }
                        },
                        "XII": {
                            "": {
                                "crystal_radius": 1.84,
                                "ionic_radius": 1.7
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "19 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 2.83,
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 0.9,
                "Youngs modulus": "no data GPa",
                "NMR Quadrupole Moment": {
                    "Ra-223": 1210.3
                },
                "Metallic radius": 2.293,
                "iupac_ordering": 12,
                "IUPAC ordering": 12
            },
            "Rb": {
                "Atomic mass": 85.4678,
                "Atomic no": 37,
                "Atomic orbitals": {
                    "1s": -540.957115,
                    "2p": -64.784678,
                    "2s": -71.291202,
                    "3d": -3.915508,
                    "3p": -8.165416,
                    "3s": -10.513861,
                    "4p": -0.59217,
                    "4s": -1.135051,
                    "5s": -0.085375
                },
                "Atomic radius": 2.35,
                "Atomic radius calculated": 2.65,
                "Boiling point": "961 K",
                "Brinell hardness": "0.216 MN m<sup>-2</sup>",
                "Bulk modulus": "2.5 GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    1
                ],
                "Critical temperature": "2093 K",
                "Density of solid": "1532 kg m<sup>-3</sup>",
                "Electrical resistivity": "13.3 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Kr].5s<sup>1</sup>",
                "ICSD oxidation states": [
                    1
                ],
                "Ionic radii": {
                    "1": 1.66
                },
                "Liquid range": "648.54 K",
                "Melting point": "312.46 K",
                "Mendeleev no": 9,
                "Mineral hardness": "0.3",
                "Molar volume": "55.76 cm<sup>3</sup>",
                "Name": "Rubidium",
                "Oxidation states": [
                    1
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "1": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.66,
                                "ionic_radius": 1.52
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.7,
                                "ionic_radius": 1.56
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.75,
                                "ionic_radius": 1.61
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.77,
                                "ionic_radius": 1.63
                            }
                        },
                        "X": {
                            "": {
                                "crystal_radius": 1.8,
                                "ionic_radius": 1.66
                            }
                        },
                        "XI": {
                            "": {
                                "crystal_radius": 1.83,
                                "ionic_radius": 1.69
                            }
                        },
                        "XII": {
                            "": {
                                "crystal_radius": 1.86,
                                "ionic_radius": 1.72
                            }
                        },
                        "XIV": {
                            "": {
                                "crystal_radius": 1.97,
                                "ionic_radius": 1.83
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "58 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 3.03,
                "Velocity of sound": "1300 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 0.82,
                "Youngs modulus": "2.4 GPa",
                "Metallic radius": 2.537,
                "iupac_ordering": 8,
                "IUPAC ordering": 8
            },
            "Re": {
                "Atomic mass": 186.207,
                "Atomic no": 75,
                "Atomic orbitals": {
                    "1s": -2407.665572,
                    "2p": -380.982869,
                    "2s": -397.087707,
                    "3d": -69.57676,
                    "3p": -83.634578,
                    "3s": -91.149193,
                    "4d": -9.516816,
                    "4f": -1.92508,
                    "4p": -15.295495,
                    "4s": -18.454325,
                    "5d": -0.258639,
                    "5p": -1.631227,
                    "5s": -2.567348,
                    "6s": -0.186859
                },
                "Atomic radius": 1.35,
                "Atomic radius calculated": 1.88,
                "Boiling point": "5869 K",
                "Brinell hardness": "1320 MN m<sup>-2</sup>",
                "Bulk modulus": "370 GPa",
                "Coefficient of linear thermal expansion": "6.2 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    4
                ],
                "Critical temperature": "no data K",
                "Density of solid": "21020 kg m<sup>-3</sup>",
                "Electrical resistivity": "18 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>14</sup>.5d<sup>5</sup>.6s<sup>2</sup>",
                "ICSD oxidation states": [
                    3,
                    4,
                    5,
                    6,
                    7
                ],
                "Ionic radii": {
                    "4": 0.77,
                    "5": 0.72,
                    "6": 0.69,
                    "7": 0.67
                },
                "Liquid range": "2410 K",
                "Melting point": "3459 K",
                "Mendeleev no": 58,
                "Mineral hardness": "7.0",
                "Molar volume": "8.86 cm<sup>3</sup>",
                "Name": "Rhenium",
                "Oxidation states": [
                    -3,
                    -1,
                    1,
                    2,
                    3,
                    4,
                    5,
                    6,
                    7
                ],
                "Poissons ratio": "0.30",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "178 GPa",
                "Shannon radii": {
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.77,
                                "ionic_radius": 0.63
                            }
                        }
                    },
                    "5": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.72,
                                "ionic_radius": 0.58
                            }
                        }
                    },
                    "6": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.69,
                                "ionic_radius": 0.55
                            }
                        }
                    },
                    "7": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.52,
                                "ionic_radius": 0.38
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.67,
                                "ionic_radius": 0.53
                            }
                        }
                    }
                },
                "Superconduction temperature": "1.70 K",
                "Thermal conductivity": "48 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "4700 m s<sup>-1</sup>",
                "Vickers hardness": "2450 MN m<sup>-2</sup>",
                "X": 1.9,
                "Youngs modulus": "463 GPa",
                "Metallic radius": 1.375,
                "iupac_ordering": 59,
                "IUPAC ordering": 59
            },
            "Rh": {
                "Atomic mass": 102.9055,
                "Atomic no": 45,
                "Atomic orbitals": {
                    "1s": -821.136773,
                    "2p": -108.357665,
                    "2s": -116.80695,
                    "3d": -11.21725,
                    "3p": -17.415299,
                    "3s": -20.765603,
                    "4d": -0.239422,
                    "4p": -1.806456,
                    "4s": -2.825505,
                    "5s": -0.154624
                },
                "Atomic radius": 1.35,
                "Atomic radius calculated": 1.73,
                "Boiling point": "3968 K",
                "Brinell hardness": "1100 MN m<sup>-2</sup>",
                "Bulk modulus": "380 GPa",
                "Coefficient of linear thermal expansion": "8.2 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "12450 kg m<sup>-3</sup>",
                "Electrical resistivity": "4.3 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Kr].4d<sup>8</sup>.5s<sup>1</sup>",
                "ICSD oxidation states": [
                    3,
                    4
                ],
                "Ionic radii": {
                    "3": 0.805,
                    "4": 0.74,
                    "5": 0.69
                },
                "Liquid range": "1731 K",
                "Melting point": "2237 K",
                "Mendeleev no": 65,
                "Mineral hardness": "6.0",
                "Molar volume": "8.28 cm<sup>3</sup>",
                "Name": "Rhodium",
                "Oxidation states": [
                    -1,
                    1,
                    2,
                    3,
                    4,
                    5,
                    6
                ],
                "Poissons ratio": "0.26",
                "Reflectivity": "84 %",
                "Refractive index": "no data",
                "Rigidity modulus": "150 GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.805,
                                "ionic_radius": 0.665
                            }
                        }
                    },
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.74,
                                "ionic_radius": 0.6
                            }
                        }
                    },
                    "5": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.69,
                                "ionic_radius": 0.55
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "150 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "4700 m s<sup>-1</sup>",
                "Vickers hardness": "1246 MN m<sup>-2</sup>",
                "X": 2.28,
                "Youngs modulus": "275 GPa",
                "Metallic radius": 1.345,
                "iupac_ordering": 66,
                "IUPAC ordering": 66
            },
            "Rn": {
                "Atomic mass": 220.0,
                "Atomic no": 86,
                "Atomic orbitals": {
                    "1s": -3204.756288,
                    "2p": -527.533025,
                    "2s": -546.57796,
                    "3d": -106.945006,
                    "3p": -124.172862,
                    "3s": -133.369144,
                    "4d": -19.449994,
                    "4f": -8.953318,
                    "4p": -27.108985,
                    "4s": -31.230804,
                    "5d": -1.911329,
                    "5p": -4.408702,
                    "5s": -5.889683,
                    "6p": -0.29318,
                    "6s": -0.62657
                },
                "Atomic radius": "no data",
                "Atomic radius calculated": 1.2,
                "Boiling point": "211.3 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Critical temperature": "377 K",
                "Density of solid": "no data kg m<sup>-3</sup>",
                "Electrical resistivity": "no data 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>14</sup>.5d<sup>10</sup>.6s<sup>2</sup>.6p<sup>6</sup>",
                "Liquid range": "9.3 K",
                "Max oxidation state": 0.0,
                "Melting point": "202 K",
                "Mendeleev no": 6,
                "Min oxidation state": 0.0,
                "Mineral hardness": "no data",
                "Molar volume": "50.50 cm<sup>3</sup>",
                "Name": "Radon",
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "0.00361 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 2.2,
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 2.2,
                "Youngs modulus": "no data GPa",
                "Metallic radius": "no data",
                "iupac_ordering": 0,
                "IUPAC ordering": 0
            },
            "Ru": {
                "Atomic mass": 101.07,
                "Atomic no": 44,
                "Atomic orbitals": {
                    "1s": -782.918621,
                    "2p": -102.333649,
                    "2s": -110.536054,
                    "3d": -10.195668,
                    "3p": -16.145217,
                    "3s": -19.366692,
                    "4d": -0.210375,
                    "4p": -1.667549,
                    "4s": -2.628363,
                    "5s": -0.152834
                },
                "Atomic radius": 1.3,
                "Atomic radius calculated": 1.78,
                "Boiling point": "4423 K",
                "Brinell hardness": "2160 MN m<sup>-2</sup>",
                "Bulk modulus": "220 GPa",
                "Coefficient of linear thermal expansion": "6.4 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3,
                    4
                ],
                "Critical temperature": "no data K",
                "Density of solid": "12370 kg m<sup>-3</sup>",
                "Electrical resistivity": "7.1 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Kr].4d<sup>7</sup>.5s<sup>1</sup>",
                "ICSD oxidation states": [
                    2,
                    3,
                    4,
                    5,
                    6
                ],
                "Ionic radii": {
                    "3": 0.82,
                    "4": 0.76,
                    "5": 0.705,
                    "7": 0.52,
                    "8": 0.5
                },
                "Liquid range": "1816 K",
                "Melting point": "2607 K",
                "Mendeleev no": 62,
                "Mineral hardness": "6.5",
                "Molar volume": "8.17 cm<sup>3</sup>",
                "Name": "Ruthenium",
                "Oxidation states": [
                    -2,
                    1,
                    2,
                    3,
                    4,
                    5,
                    6,
                    7,
                    8
                ],
                "Poissons ratio": "0.30",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "173 GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.82,
                                "ionic_radius": 0.68
                            }
                        }
                    },
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.76,
                                "ionic_radius": 0.62
                            }
                        }
                    },
                    "5": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.705,
                                "ionic_radius": 0.565
                            }
                        }
                    },
                    "7": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.52,
                                "ionic_radius": 0.38
                            }
                        }
                    },
                    "8": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.5,
                                "ionic_radius": 0.36
                            }
                        }
                    }
                },
                "Superconduction temperature": "0.49 K",
                "Thermal conductivity": "120 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "5970 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 2.2,
                "Youngs modulus": "447 GPa",
                "Metallic radius": 1.339,
                "iupac_ordering": 63,
                "IUPAC ordering": 63
            },
            "S": {
                "Atomic mass": 32.065,
                "Atomic no": 16,
                "Atomic orbitals": {
                    "1s": -87.789937,
                    "2p": -5.751257,
                    "2s": -7.69994,
                    "3p": -0.261676,
                    "3s": -0.630912
                },
                "Atomic radius": 1.0,
                "Atomic radius calculated": 0.88,
                "Boiling point": "717.87 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "7.7 GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    -2,
                    2,
                    4,
                    6
                ],
                "Critical temperature": "1314 K",
                "Density of solid": "1960 kg m<sup>-3</sup>",
                "Electrical resistivity": "&gt; 10<sup>23</sup>10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ne].3s<sup>2</sup>.3p<sup>4</sup>",
                "ICSD oxidation states": [
                    -1,
                    2,
                    4,
                    -2,
                    6
                ],
                "Ionic radii": {
                    "-2": 1.7,
                    "4": 0.51,
                    "6": 0.43
                },
                "Liquid range": "329.51 K",
                "Melting point": "388.36 K",
                "Mendeleev no": 94,
                "Mineral hardness": "2.0",
                "Molar volume": "15.53 cm<sup>3</sup>",
                "Name": "Sulfur",
                "Oxidation states": [
                    -2,
                    -1,
                    1,
                    2,
                    3,
                    4,
                    5,
                    6
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "1.001111",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "-2": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.7,
                                "ionic_radius": 1.84
                            }
                        }
                    },
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.51,
                                "ionic_radius": 0.37
                            }
                        }
                    },
                    "6": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.26,
                                "ionic_radius": 0.12
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.43,
                                "ionic_radius": 0.29
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "0.205 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.8,
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 2.58,
                "Youngs modulus": "no data GPa",
                "NMR Quadrupole Moment": {
                    "S-33": -67.8,
                    "S-35": 47.1
                },
                "Metallic radius": "no data",
                "iupac_ordering": 96,
                "IUPAC ordering": 96
            },
            "Sb": {
                "Atomic mass": 121.76,
                "Atomic no": 51,
                "Atomic orbitals": {
                    "1s": -1070.823495,
                    "2p": -149.214271,
                    "2s": -159.171745,
                    "3d": -19.239895,
                    "3p": -26.956184,
                    "3s": -31.098242,
                    "4d": -1.297338,
                    "4p": -3.646579,
                    "4s": -5.04964,
                    "5p": -0.185623,
                    "5s": -0.445605
                },
                "Atomic radius": 1.45,
                "Atomic radius calculated": 1.33,
                "Boiling point": "1860 K",
                "Brinell hardness": "294 MN m<sup>-2</sup>",
                "Bulk modulus": "42 GPa",
                "Coefficient of linear thermal expansion": "11 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    -3,
                    3,
                    5
                ],
                "Critical temperature": "no data K",
                "Density of solid": "6697 kg m<sup>-3</sup>",
                "Electrical resistivity": "40 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Kr].4d<sup>10</sup>.5s<sup>2</sup>.5p<sup>3</sup>",
                "ICSD oxidation states": [
                    -2,
                    3,
                    5,
                    -1,
                    -3
                ],
                "Ionic radii": {
                    "3": 0.9,
                    "5": 0.76
                },
                "Liquid range": "956.22 K",
                "Melting point": "903.78 K",
                "Mendeleev no": 88,
                "Mineral hardness": "3.0",
                "Molar volume": "18.19 cm<sup>3</sup>",
                "Name": "Antimony",
                "Oxidation states": [
                    -3,
                    3,
                    5
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "55 %",
                "Refractive index": "no data",
                "Rigidity modulus": "20 GPa",
                "Shannon radii": {
                    "3": {
                        "IVPY": {
                            "": {
                                "crystal_radius": 0.9,
                                "ionic_radius": 0.76
                            }
                        },
                        "V": {
                            "": {
                                "crystal_radius": 0.94,
                                "ionic_radius": 0.8
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.9,
                                "ionic_radius": 0.76
                            }
                        }
                    },
                    "5": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.74,
                                "ionic_radius": 0.6
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "24 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 2.06,
                "Velocity of sound": "3420 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 2.05,
                "Youngs modulus": "55 GPa",
                "NMR Quadrupole Moment": {
                    "Sb-121": -543.11,
                    "Sb-123": -692.14
                },
                "Metallic radius": 1.61,
                "iupac_ordering": 88,
                "IUPAC ordering": 88
            },
            "Sc": {
                "Atomic mass": 44.955912,
                "Atomic no": 21,
                "Atomic orbitals": {
                    "1s": -160.184109,
                    "2p": -14.240006,
                    "2s": -17.206464,
                    "3d": -0.13108,
                    "3p": -1.233165,
                    "3s": -1.988378,
                    "4s": -0.156478
                },
                "Atomic radius": 1.6,
                "Atomic radius calculated": 1.84,
                "Boiling point": "3103 K",
                "Brinell hardness": "750 MN m<sup>-2</sup>",
                "Bulk modulus": "57 GPa",
                "Coefficient of linear thermal expansion": "10.2 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "2985 kg m<sup>-3</sup>",
                "Electrical resistivity": "about 55 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ar].3d<sup>1</sup>.4s<sup>2</sup>",
                "ICSD oxidation states": [
                    2,
                    3
                ],
                "Ionic radii": {
                    "3": 0.885
                },
                "Liquid range": "1289 K",
                "Melting point": "1814 K",
                "Mendeleev no": 19,
                "Mineral hardness": "no data",
                "Molar volume": "15.00 cm<sup>3</sup>",
                "Name": "Scandium",
                "Oxidation states": [
                    1,
                    2,
                    3
                ],
                "Poissons ratio": "0.28",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "29 GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.885,
                                "ionic_radius": 0.745
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.01,
                                "ionic_radius": 0.87
                            }
                        }
                    }
                },
                "Superconduction temperature": "0.05  (under pressure)K",
                "Thermal conductivity": "16 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 2.11,
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.36,
                "Youngs modulus": "74 GPa",
                "NMR Quadrupole Moment": {
                    "Sc-45": -220.2
                },
                "Metallic radius": 1.641,
                "iupac_ordering": 49,
                "IUPAC ordering": 49
            },
            "Se": {
                "Atomic mass": 78.96,
                "Atomic no": 34,
                "Atomic orbitals": {
                    "1s": -451.300258,
                    "2p": -51.514388,
                    "2s": -57.311948,
                    "3d": -2.011392,
                    "3p": -5.553517,
                    "3s": -7.547186,
                    "4p": -0.245806,
                    "4s": -0.621248
                },
                "Atomic radius": 1.15,
                "Atomic radius calculated": 1.03,
                "Boiling point": "958 K",
                "Brinell hardness": "736 MN m<sup>-2</sup>",
                "Bulk modulus": "8.3 GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    -2,
                    2,
                    4,
                    6
                ],
                "Critical temperature": "1766 K",
                "Density of solid": "4819 kg m<sup>-3</sup>",
                "Electrical resistivity": "high 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ar].3d<sup>10</sup>.4s<sup>2</sup>.4p<sup>4</sup>",
                "ICSD oxidation states": [
                    -1,
                    4,
                    -2,
                    6
                ],
                "Ionic radii": {
                    "-2": 1.84,
                    "4": 0.64,
                    "6": 0.56
                },
                "Liquid range": "464 K",
                "Melting point": "494 K",
                "Mendeleev no": 93,
                "Mineral hardness": "2.0",
                "Molar volume": "16.42 cm<sup>3</sup>",
                "Name": "Selenium",
                "Oxidation states": [
                    -2,
                    2,
                    4,
                    6
                ],
                "Poissons ratio": "0.33",
                "Reflectivity": "no data %",
                "Refractive index": "1.000895",
                "Rigidity modulus": "3.7 GPa",
                "Shannon radii": {
                    "-2": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.84,
                                "ionic_radius": 1.98
                            }
                        }
                    },
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.64,
                                "ionic_radius": 0.5
                            }
                        }
                    },
                    "6": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.42,
                                "ionic_radius": 0.28
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.56,
                                "ionic_radius": 0.42
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "0.52 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.9,
                "Velocity of sound": "3350 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 2.55,
                "Youngs modulus": "10 GPa",
                "Metallic radius": "no data",
                "iupac_ordering": 95,
                "IUPAC ordering": 95
            },
            "Si": {
                "Atomic mass": 28.0855,
                "Atomic no": 14,
                "Atomic orbitals": {
                    "1s": -65.184426,
                    "2p": -3.514938,
                    "2s": -5.075056,
                    "3p": -0.153293,
                    "3s": -0.398139
                },
                "Atomic radius": 1.1,
                "Atomic radius calculated": 1.11,
                "Boiling point": "3173 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "100 GPa",
                "Coefficient of linear thermal expansion": "2.6 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    -4,
                    4
                ],
                "Critical temperature": "no data K",
                "Density of solid": "2330 kg m<sup>-3</sup>",
                "Electrical resistivity": "about 100000 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ne].3s<sup>2</sup>.3p<sup>2</sup>",
                "ICSD oxidation states": [
                    -4,
                    4
                ],
                "Ionic radii": {
                    "4": 0.54
                },
                "Liquid range": "1486 K",
                "Melting point": "1687 K",
                "Mendeleev no": 85,
                "Mineral hardness": "6.5",
                "Molar volume": "12.06 cm<sup>3</sup>",
                "Name": "Silicon",
                "Oxidation states": [
                    -4,
                    -3,
                    -2,
                    -1,
                    1,
                    2,
                    3,
                    4
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "28 %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "4": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.4,
                                "ionic_radius": 0.26
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.54,
                                "ionic_radius": 0.4
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "150 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 2.1,
                "Velocity of sound": "2200 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.9,
                "Youngs modulus": "47 GPa",
                "Metallic radius": "no data",
                "iupac_ordering": 85,
                "IUPAC ordering": 85
            },
            "Sm": {
                "Atomic mass": 150.36,
                "Atomic no": 62,
                "Atomic orbitals": {
                    "1s": -1617.183426,
                    "2p": -242.729726,
                    "2s": -255.498846,
                    "3d": -39.528656,
                    "3p": -50.08426,
                    "3s": -55.731133,
                    "4d": -4.814978,
                    "4f": -0.21776,
                    "4p": -8.672685,
                    "4s": -10.844667,
                    "5p": -0.835987,
                    "5s": -1.408552,
                    "6s": -0.128259
                },
                "Atomic radius": 1.85,
                "Atomic radius calculated": 2.38,
                "Boiling point": "2076 K",
                "Brinell hardness": "441 MN m<sup>-2</sup>",
                "Bulk modulus": "38 GPa",
                "Coefficient of linear thermal expansion": "12.7 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "7353 kg m<sup>-3</sup>",
                "Electrical resistivity": "94 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>6</sup>.6s<sup>2</sup>",
                "ICSD oxidation states": [
                    2,
                    3
                ],
                "Ionic radii": {
                    "2": 1.36,
                    "3": 1.0979999999999999
                },
                "Liquid range": "731 K",
                "Melting point": "1345 K",
                "Mendeleev no": 28,
                "Mineral hardness": "no data",
                "Molar volume": "19.98 cm<sup>3</sup>",
                "Name": "Samarium",
                "Oxidation states": [
                    2,
                    3
                ],
                "Poissons ratio": "0.27",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "20 GPa",
                "Shannon radii": {
                    "2": {
                        "VII": {
                            "": {
                                "crystal_radius": 1.36,
                                "ionic_radius": 1.22
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.41,
                                "ionic_radius": 1.27
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.46,
                                "ionic_radius": 1.32
                            }
                        }
                    },
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.098,
                                "ionic_radius": 0.958
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.16,
                                "ionic_radius": 1.02
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.219,
                                "ionic_radius": 1.079
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.272,
                                "ionic_radius": 1.132
                            }
                        },
                        "XII": {
                            "": {
                                "crystal_radius": 1.38,
                                "ionic_radius": 1.24
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "13 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "2130 m s<sup>-1</sup>",
                "Vickers hardness": "412 MN m<sup>-2</sup>",
                "X": 1.17,
                "Youngs modulus": "50 GPa",
                "Metallic radius": 1.804,
                "iupac_ordering": 42,
                "IUPAC ordering": 42
            },
            "Sn": {
                "Atomic mass": 118.71,
                "Atomic no": 50,
                "Atomic orbitals": {
                    "1s": -1026.762169,
                    "2p": -141.821093,
                    "2s": -151.523991,
                    "3d": -17.657276,
                    "3p": -25.117913,
                    "3s": -29.125969,
                    "4d": -1.004952,
                    "4p": -3.211998,
                    "4s": -4.546335,
                    "5p": -0.14445,
                    "5s": -0.369349
                },
                "Atomic radius": 1.45,
                "Atomic radius calculated": 1.45,
                "Boiling point": "2875 K",
                "Brinell hardness": "51 MN m<sup>-2</sup>",
                "Bulk modulus": "58 GPa",
                "Coefficient of linear thermal expansion": "22 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    -4,
                    2,
                    4
                ],
                "Critical temperature": "no data K",
                "Density of solid": "7310 kg m<sup>-3</sup>",
                "Electrical resistivity": "11.5 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Kr].4d<sup>10</sup>.5s<sup>2</sup>.5p<sup>2</sup>",
                "ICSD oxidation states": [
                    2,
                    3,
                    4
                ],
                "Ionic radii": {
                    "4": 0.83
                },
                "Liquid range": "2369.92 K",
                "Melting point": "505.08 K",
                "Mendeleev no": 83,
                "Mineral hardness": "1.5",
                "Molar volume": "16.29 cm<sup>3</sup>",
                "Name": "Tin",
                "Oxidation states": [
                    -4,
                    2,
                    4
                ],
                "Poissons ratio": "0.36",
                "Reflectivity": "54 %",
                "Refractive index": "no data",
                "Rigidity modulus": "18 GPa",
                "Shannon radii": {
                    "4": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.69,
                                "ionic_radius": 0.55
                            }
                        },
                        "V": {
                            "": {
                                "crystal_radius": 0.76,
                                "ionic_radius": 0.62
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.83,
                                "ionic_radius": 0.69
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 0.89,
                                "ionic_radius": 0.75
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 0.95,
                                "ionic_radius": 0.81
                            }
                        }
                    }
                },
                "Superconduction temperature": "3.72 K",
                "Thermal conductivity": "67 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 2.17,
                "Velocity of sound": "2500 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.96,
                "Youngs modulus": "50 GPa",
                "NMR Quadrupole Moment": {
                    "Sn-119": -132.1
                },
                "Metallic radius": 1.58,
                "iupac_ordering": 83,
                "IUPAC ordering": 83
            },
            "Sr": {
                "Atomic mass": 87.62,
                "Atomic no": 38,
                "Atomic orbitals": {
                    "1s": -572.870169,
                    "2p": -69.745941,
                    "2s": -76.491823,
                    "3d": -4.813498,
                    "3p": -9.301863,
                    "3s": -11.771585,
                    "4p": -0.844489,
                    "4s": -1.455317,
                    "5s": -0.131793
                },
                "Atomic radius": 2.0,
                "Atomic radius calculated": 2.19,
                "Boiling point": "1655 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "22.5 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    2
                ],
                "Critical temperature": "no data K",
                "Density of solid": "2630 kg m<sup>-3</sup>",
                "Electrical resistivity": "13.5 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Kr].5s<sup>2</sup>",
                "ICSD oxidation states": [
                    2
                ],
                "Ionic radii": {
                    "2": 1.32
                },
                "Liquid range": "605 K",
                "Melting point": "1050 K",
                "Mendeleev no": 15,
                "Mineral hardness": "1.5",
                "Molar volume": "33.94 cm<sup>3</sup>",
                "Name": "Strontium",
                "Oxidation states": [
                    2
                ],
                "Poissons ratio": "0.28",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "6.1 GPa",
                "Shannon radii": {
                    "2": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.32,
                                "ionic_radius": 1.18
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.35,
                                "ionic_radius": 1.21
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.4,
                                "ionic_radius": 1.26
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.45,
                                "ionic_radius": 1.31
                            }
                        },
                        "X": {
                            "": {
                                "crystal_radius": 1.5,
                                "ionic_radius": 1.36
                            }
                        },
                        "XII": {
                            "": {
                                "crystal_radius": 1.58,
                                "ionic_radius": 1.44
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "35 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 2.49,
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 0.95,
                "Youngs modulus": "no data GPa",
                "NMR Quadrupole Moment": {
                    "Sr-87": 305.2
                },
                "Metallic radius": 2.151,
                "iupac_ordering": 14,
                "IUPAC ordering": 14
            },
            "Ta": {
                "Atomic mass": 180.94788,
                "Atomic no": 73,
                "Atomic orbitals": {
                    "1s": -2275.371387,
                    "2p": -357.248334,
                    "2s": -372.828724,
                    "3d": -63.942521,
                    "3p": -77.440942,
                    "3s": -84.658467,
                    "4d": -8.265848,
                    "4f": -1.199347,
                    "4p": -13.71981,
                    "4s": -16.713337,
                    "5d": -0.182464,
                    "5p": -1.37653,
                    "5s": -2.223807,
                    "6s": -0.174814
                },
                "Atomic radius": 1.45,
                "Atomic radius calculated": 2.0,
                "Boiling point": "5731 K",
                "Brinell hardness": "800 MN m<sup>-2</sup>",
                "Bulk modulus": "200 GPa",
                "Coefficient of linear thermal expansion": "6.3 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    5
                ],
                "Critical temperature": "no data K",
                "Density of solid": "16650 kg m<sup>-3</sup>",
                "Electrical resistivity": "13.5 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>14</sup>.5d<sup>3</sup>.6s<sup>2</sup>",
                "ICSD oxidation states": [
                    3,
                    4,
                    5
                ],
                "Ionic radii": {
                    "3": 0.86,
                    "4": 0.82,
                    "5": 0.78
                },
                "Liquid range": "2441 K",
                "Melting point": "3290 K",
                "Mendeleev no": 52,
                "Mineral hardness": "6.5",
                "Molar volume": "10.85 cm<sup>3</sup>",
                "Name": "Tantalum",
                "Oxidation states": [
                    -1,
                    2,
                    3,
                    4,
                    5
                ],
                "Poissons ratio": "0.34",
                "Reflectivity": "78 %",
                "Refractive index": "no data",
                "Rigidity modulus": "69 GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.86,
                                "ionic_radius": 0.72
                            }
                        }
                    },
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.82,
                                "ionic_radius": 0.68
                            }
                        }
                    },
                    "5": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.78,
                                "ionic_radius": 0.64
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 0.83,
                                "ionic_radius": 0.69
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 0.88,
                                "ionic_radius": 0.74
                            }
                        }
                    }
                },
                "Superconduction temperature": "4.47 K",
                "Thermal conductivity": "57 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "3400 m s<sup>-1</sup>",
                "Vickers hardness": "873 MN m<sup>-2</sup>",
                "X": 1.5,
                "Youngs modulus": "186 GPa",
                "Metallic radius": 1.47,
                "iupac_ordering": 53,
                "IUPAC ordering": 53
            },
            "Tb": {
                "Atomic mass": 158.92535,
                "Atomic no": 65,
                "Atomic orbitals": {
                    "1s": -1785.331942,
                    "2p": -271.590585,
                    "2s": -285.121013,
                    "3d": -45.443863,
                    "3p": -56.785113,
                    "3s": -62.851563,
                    "4d": -5.467662,
                    "4f": -0.256311,
                    "4p": -9.735637,
                    "4s": -12.120486,
                    "5p": -0.88723,
                    "5s": -1.513669,
                    "6s": -0.131677
                },
                "Atomic radius": 1.75,
                "Atomic radius calculated": 2.25,
                "Boiling point": "3503 K",
                "Brinell hardness": "677 MN m<sup>-2</sup>",
                "Bulk modulus": "38.7 GPa",
                "Coefficient of linear thermal expansion": "10.3 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "8219 kg m<sup>-3</sup>",
                "Electrical resistivity": "115 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>9</sup>.6s<sup>2</sup>",
                "ICSD oxidation states": [
                    3,
                    4
                ],
                "Ionic radii": {
                    "3": 1.063,
                    "4": 0.9
                },
                "Liquid range": "1874 K",
                "Melting point": "1629 K",
                "Mendeleev no": 26,
                "Mineral hardness": "no data",
                "Molar volume": "19.30 cm<sup>3</sup>",
                "Name": "Terbium",
                "Oxidation states": [
                    1,
                    3,
                    4
                ],
                "Poissons ratio": "0.26",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "22 GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.063,
                                "ionic_radius": 0.923
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.12,
                                "ionic_radius": 0.98
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.18,
                                "ionic_radius": 1.04
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.235,
                                "ionic_radius": 1.095
                            }
                        }
                    },
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.9,
                                "ionic_radius": 0.76
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.02,
                                "ionic_radius": 0.88
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "11 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "2620 m s<sup>-1</sup>",
                "Vickers hardness": "863 MN m<sup>-2</sup>",
                "X": 1.1,
                "Youngs modulus": "56 GPa",
                "Metallic radius": 1.781,
                "iupac_ordering": 39,
                "IUPAC ordering": 39
            },
            "Tc": {
                "Atomic mass": 98.0,
                "Atomic no": 43,
                "Atomic orbitals": {
                    "1s": -745.742024,
                    "2p": -96.61021,
                    "2s": -104.567508,
                    "3d": -9.33986,
                    "3p": -15.041738,
                    "3s": -18.135303,
                    "4d": -0.270262,
                    "4p": -1.64323,
                    "4s": -2.550712,
                    "5s": -0.183636
                },
                "Atomic radius": 1.35,
                "Atomic radius calculated": 1.83,
                "Boiling point": "4538 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    4,
                    7
                ],
                "Critical temperature": "no data K",
                "Density of solid": "11500 kg m<sup>-3</sup>",
                "Electrical resistivity": "about 22 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Kr].4d<sup>5</sup>.5s<sup>2</sup>",
                "Ionic radii": {
                    "4": 0.785,
                    "5": 0.74,
                    "7": 0.7
                },
                "Liquid range": "2108 K",
                "Melting point": "2430 K",
                "Mendeleev no": 59,
                "Mineral hardness": "no data",
                "Molar volume": "8.63 cm<sup>3</sup>",
                "Name": "Technetium",
                "Oxidation states": [
                    -3,
                    -1,
                    1,
                    2,
                    3,
                    4,
                    5,
                    6,
                    7
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.785,
                                "ionic_radius": 0.645
                            }
                        }
                    },
                    "5": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.74,
                                "ionic_radius": 0.6
                            }
                        }
                    },
                    "7": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.51,
                                "ionic_radius": 0.37
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.7,
                                "ionic_radius": 0.56
                            }
                        }
                    }
                },
                "Superconduction temperature": "7.8 K",
                "Thermal conductivity": "51 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.9,
                "Youngs modulus": "no data GPa",
                "Metallic radius": 1.363,
                "iupac_ordering": 60,
                "IUPAC ordering": 60
            },
            "Te": {
                "Atomic mass": 127.6,
                "Atomic no": 52,
                "Atomic orbitals": {
                    "1s": -1115.831819,
                    "2p": -156.808583,
                    "2s": -167.021776,
                    "3d": -20.887801,
                    "3p": -28.860685,
                    "3s": -33.137485,
                    "4d": -1.608381,
                    "4p": -4.100084,
                    "4s": -5.572846,
                    "5p": -0.226594,
                    "5s": -0.520997
                },
                "Atomic radius": 1.4,
                "Atomic radius calculated": 1.23,
                "Boiling point": "1261 K",
                "Brinell hardness": "180 MN m<sup>-2</sup>",
                "Bulk modulus": "65 GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    -2,
                    2,
                    4,
                    6
                ],
                "Critical temperature": "no data K",
                "Density of solid": "6240 kg m<sup>-3</sup>",
                "Electrical resistivity": "about 10000 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Kr].4d<sup>10</sup>.5s<sup>2</sup>.5p<sup>4</sup>",
                "ICSD oxidation states": [
                    -2,
                    4,
                    -1,
                    6
                ],
                "Ionic radii": {
                    "-2": 2.07,
                    "4": 1.11,
                    "6": 0.7
                },
                "Liquid range": "538.34 K",
                "Melting point": "722.66 K",
                "Mendeleev no": 92,
                "Mineral hardness": "2.25",
                "Molar volume": "20.46 cm<sup>3</sup>",
                "Name": "Tellurium",
                "Oxidation states": [
                    -2,
                    2,
                    4,
                    5,
                    6
                ],
                "Poissons ratio": "no data",
                "Reflectivity": "50 %",
                "Refractive index": "1.000991",
                "Rigidity modulus": "16 GPa",
                "Shannon radii": {
                    "-2": {
                        "VI": {
                            "": {
                                "crystal_radius": 2.07,
                                "ionic_radius": 2.21
                            }
                        }
                    },
                    "4": {
                        "III": {
                            "": {
                                "crystal_radius": 0.66,
                                "ionic_radius": 0.52
                            }
                        },
                        "IV": {
                            "": {
                                "crystal_radius": 0.8,
                                "ionic_radius": 0.66
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 1.11,
                                "ionic_radius": 0.97
                            }
                        }
                    },
                    "6": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.57,
                                "ionic_radius": 0.43
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.7,
                                "ionic_radius": 0.56
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "3 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 2.06,
                "Velocity of sound": "2610 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 2.1,
                "Youngs modulus": "43 GPa",
                "Metallic radius": "no data",
                "iupac_ordering": 94,
                "IUPAC ordering": 94
            },
            "Th": {
                "Atomic mass": 232.03806,
                "Atomic no": 90,
                "Atomic orbitals": {
                    "1s": -3524.439052,
                    "2p": -588.218112,
                    "2s": -608.350958,
                    "3d": -123.846396,
                    "3p": -142.25581,
                    "3s": -152.079741,
                    "4d": -24.955184,
                    "4f": -13.397389,
                    "4p": -33.325252,
                    "4s": -37.814094,
                    "5d": -3.625729,
                    "5p": -6.58281,
                    "5s": -8.287057,
                    "6d": -0.172896,
                    "6p": -0.846921,
                    "6s": -1.333769,
                    "7s": -0.135872
                },
                "Atomic radius": 1.8,
                "Atomic radius calculated": "no data",
                "Boiling point": "5093 K",
                "Brinell hardness": "400 MN m<sup>-2</sup>",
                "Bulk modulus": "54 GPa",
                "Coefficient of linear thermal expansion": "11.0 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    4
                ],
                "Critical temperature": "no data K",
                "Density of solid": "11724 kg m<sup>-3</sup>",
                "Electrical resistivity": "15 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Rn].6d<sup>2</sup>.7s<sup>2</sup>",
                "ICSD oxidation states": [
                    4
                ],
                "Ionic radii": {
                    "4": 1.08
                },
                "Liquid range": "2978 K",
                "Melting point": "2115 K",
                "Mendeleev no": 47,
                "Mineral hardness": "3.0",
                "Molar volume": "19.80 cm<sup>3</sup>",
                "Name": "Thorium",
                "Oxidation states": [
                    2,
                    3,
                    4
                ],
                "Poissons ratio": "0.27",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "31 GPa",
                "Shannon radii": {
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.08,
                                "ionic_radius": 0.94
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.19,
                                "ionic_radius": 1.05
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.23,
                                "ionic_radius": 1.09
                            }
                        },
                        "X": {
                            "": {
                                "crystal_radius": 1.27,
                                "ionic_radius": 1.13
                            }
                        },
                        "XI": {
                            "": {
                                "crystal_radius": 1.32,
                                "ionic_radius": 1.18
                            }
                        },
                        "XII": {
                            "": {
                                "crystal_radius": 1.35,
                                "ionic_radius": 1.21
                            }
                        }
                    }
                },
                "Superconduction temperature": "1.38 K",
                "Thermal conductivity": "54 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "2490 m s<sup>-1</sup>",
                "Vickers hardness": "350 MN m<sup>-2</sup>",
                "X": 1.3,
                "Youngs modulus": "79 GPa",
                "Metallic radius": 1.798,
                "iupac_ordering": 31,
                "IUPAC ordering": 31
            },
            "Ti": {
                "Atomic mass": 47.867,
                "Atomic no": 22,
                "Atomic orbitals": {
                    "1s": -177.276643,
                    "2p": -16.285339,
                    "2s": -19.457901,
                    "3d": -0.17001,
                    "3p": -1.422947,
                    "3s": -2.258007,
                    "4s": -0.167106
                },
                "Atomic radius": 1.4,
                "Atomic radius calculated": 1.76,
                "Boiling point": "3560 K",
                "Brinell hardness": "716 MN m<sup>-2</sup>",
                "Bulk modulus": "110 GPa",
                "Coefficient of linear thermal expansion": "8.6 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    4
                ],
                "Critical temperature": "no data K",
                "Density of solid": "4507 kg m<sup>-3</sup>",
                "Electrical resistivity": "about 40 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ar].3d<sup>2</sup>.4s<sup>2</sup>",
                "ICSD oxidation states": [
                    2,
                    3,
                    4
                ],
                "Ionic radii": {
                    "2": 1.0,
                    "3": 0.81,
                    "4": 0.745
                },
                "Liquid range": "1619 K",
                "Melting point": "1941 K",
                "Mendeleev no": 51,
                "Mineral hardness": "6.0",
                "Molar volume": "10.64 cm<sup>3</sup>",
                "Name": "Titanium",
                "Oxidation states": [
                    -1,
                    2,
                    3,
                    4
                ],
                "Poissons ratio": "0.32",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "44 GPa",
                "Shannon radii": {
                    "2": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.0,
                                "ionic_radius": 0.86
                            }
                        }
                    },
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.81,
                                "ionic_radius": 0.67
                            }
                        }
                    },
                    "4": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.56,
                                "ionic_radius": 0.42
                            }
                        },
                        "V": {
                            "": {
                                "crystal_radius": 0.65,
                                "ionic_radius": 0.51
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.745,
                                "ionic_radius": 0.605
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 0.88,
                                "ionic_radius": 0.74
                            }
                        }
                    }
                },
                "Superconduction temperature": "0.40 K",
                "Thermal conductivity": "22 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "4140 m s<sup>-1</sup>",
                "Vickers hardness": "970 MN m<sup>-2</sup>",
                "X": 1.54,
                "Youngs modulus": "116 GPa",
                "NMR Quadrupole Moment": {
                    "Ti-47": 302.1,
                    "Ti-49": 247.11
                },
                "Metallic radius": 1.462,
                "iupac_ordering": 52,
                "IUPAC ordering": 52
            },
            "Tl": {
                "Atomic mass": 204.3833,
                "Atomic no": 81,
                "Atomic orbitals": {
                    "1s": -2827.569408,
                    "2p": -457.255971,
                    "2s": -474.953368,
                    "3d": -88.328299,
                    "3p": -104.099296,
                    "3s": -112.52218,
                    "4d": -14.008848,
                    "4f": -4.835747,
                    "4p": -20.797078,
                    "4s": -24.471512,
                    "5d": -0.674544,
                    "5p": -2.59873,
                    "5s": -3.811512,
                    "6p": -0.101507,
                    "6s": -0.28502
                },
                "Atomic radius": 1.9,
                "Atomic radius calculated": 1.56,
                "Boiling point": "1746 K",
                "Brinell hardness": "26.4 MN m<sup>-2</sup>",
                "Bulk modulus": "43 GPa",
                "Coefficient of linear thermal expansion": "29.9 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    1,
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "11850 kg m<sup>-3</sup>",
                "Electrical resistivity": "15 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>14</sup>.5d<sup>10</sup>.6s<sup>2</sup>.6p<sup>1</sup>",
                "ICSD oxidation states": [
                    1,
                    3
                ],
                "Ionic radii": {
                    "1": 1.64,
                    "3": 1.025
                },
                "Liquid range": "1169 K",
                "Melting point": "577 K",
                "Mendeleev no": 78,
                "Mineral hardness": "1.2",
                "Molar volume": "17.22 cm<sup>3</sup>",
                "Name": "Thallium",
                "Oxidation states": [
                    1,
                    3
                ],
                "Poissons ratio": "0.45",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "2.8 GPa",
                "Shannon radii": {
                    "1": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.64,
                                "ionic_radius": 1.5
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.73,
                                "ionic_radius": 1.59
                            }
                        },
                        "XII": {
                            "": {
                                "crystal_radius": 1.84,
                                "ionic_radius": 1.7
                            }
                        }
                    },
                    "3": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.89,
                                "ionic_radius": 0.75
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 1.025,
                                "ionic_radius": 0.885
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.12,
                                "ionic_radius": 0.98
                            }
                        }
                    }
                },
                "Superconduction temperature": "2.38 K",
                "Thermal conductivity": "46 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.96,
                "Velocity of sound": "818 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.62,
                "Youngs modulus": "8 GPa",
                "Metallic radius": 1.7,
                "iupac_ordering": 77,
                "IUPAC ordering": 77
            },
            "Tm": {
                "Atomic mass": 168.93421,
                "Atomic no": 69,
                "Atomic orbitals": {
                    "1s": -2022.471608,
                    "2p": -312.510608,
                    "2s": -327.05712,
                    "3d": -53.835494,
                    "3p": -66.239338,
                    "3s": -72.873753,
                    "4d": -6.350307,
                    "4f": -0.28312,
                    "4p": -11.187151,
                    "4s": -13.865665,
                    "5p": -0.950748,
                    "5s": -1.64999,
                    "6s": -0.135953
                },
                "Atomic radius": 1.75,
                "Atomic radius calculated": 2.22,
                "Boiling point": "2223 K",
                "Brinell hardness": "471 MN m<sup>-2</sup>",
                "Bulk modulus": "45 GPa",
                "Coefficient of linear thermal expansion": "13.3 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "9321 kg m<sup>-3</sup>",
                "Electrical resistivity": "67.6 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>13</sup>.6s<sup>2</sup>",
                "ICSD oxidation states": [
                    3
                ],
                "Ionic radii": {
                    "2": 1.17,
                    "3": 1.02
                },
                "Liquid range": "405 K",
                "Melting point": "1818 K",
                "Mendeleev no": 21,
                "Mineral hardness": "no data",
                "Molar volume": "19.1 cm<sup>3</sup>",
                "Name": "Thulium",
                "Oxidation states": [
                    2,
                    3
                ],
                "Poissons ratio": "0.21",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "31 GPa",
                "Shannon radii": {
                    "2": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.17,
                                "ionic_radius": 1.03
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.23,
                                "ionic_radius": 1.09
                            }
                        }
                    },
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.02,
                                "ionic_radius": 0.88
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.134,
                                "ionic_radius": 0.994
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.192,
                                "ionic_radius": 1.052
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "17 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "no data m s<sup>-1</sup>",
                "Vickers hardness": "520 MN m<sup>-2</sup>",
                "X": 1.25,
                "Youngs modulus": "74 GPa",
                "Metallic radius": 1.747,
                "iupac_ordering": 35,
                "IUPAC ordering": 35
            },
            "U": {
                "Atomic mass": 238.02891,
                "Atomic no": 92,
                "Atomic orbitals": {
                    "1s": -3689.355141,
                    "2p": -619.10855,
                    "2s": -639.778728,
                    "3d": -131.977358,
                    "3p": -150.97898,
                    "3s": -161.118073,
                    "4d": -27.123212,
                    "4f": -15.02746,
                    "4p": -35.853321,
                    "4s": -40.528084,
                    "5d": -3.866175,
                    "5f": -0.366543,
                    "5p": -7.018092,
                    "5s": -8.824089,
                    "6d": -0.14319,
                    "6p": -0.822538,
                    "6s": -1.325976,
                    "7s": -0.130948
                },
                "Atomic radius": 1.75,
                "Atomic radius calculated": "no data",
                "Boiling point": "4200 K",
                "Brinell hardness": "2400 MN m<sup>-2</sup>",
                "Bulk modulus": "100 GPa",
                "Coefficient of linear thermal expansion": "13.9 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    6
                ],
                "Critical temperature": "no data K",
                "Density of solid": "19050 kg m<sup>-3</sup>",
                "Electrical resistivity": "28 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Rn].5f<sup>3</sup>.6d<sup>1</sup>.7s<sup>2</sup>",
                "ICSD oxidation states": [
                    3,
                    4,
                    5,
                    6
                ],
                "Ionic radii": {
                    "3": 1.165,
                    "4": 1.03,
                    "5": 0.9,
                    "6": 0.87
                },
                "Liquid range": "2794.7 K",
                "Melting point": "1405.3 K",
                "Mendeleev no": 45,
                "Mineral hardness": "6.0",
                "Molar volume": "12.49 cm<sup>3</sup>",
                "Name": "Uranium",
                "Oxidation states": [
                    3,
                    4,
                    5,
                    6
                ],
                "Poissons ratio": "0.23",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "111 GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.165,
                                "ionic_radius": 1.025
                            }
                        }
                    },
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.03,
                                "ionic_radius": 0.89
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.09,
                                "ionic_radius": 0.95
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.14,
                                "ionic_radius": 1.0
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.19,
                                "ionic_radius": 1.05
                            }
                        },
                        "XII": {
                            "": {
                                "crystal_radius": 1.31,
                                "ionic_radius": 1.17
                            }
                        }
                    },
                    "5": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.9,
                                "ionic_radius": 0.76
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 0.98,
                                "ionic_radius": 0.84
                            }
                        }
                    },
                    "6": {
                        "II": {
                            "": {
                                "crystal_radius": 0.59,
                                "ionic_radius": 0.45
                            }
                        },
                        "IV": {
                            "": {
                                "crystal_radius": 0.66,
                                "ionic_radius": 0.52
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.87,
                                "ionic_radius": 0.73
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 0.95,
                                "ionic_radius": 0.81
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.0,
                                "ionic_radius": 0.86
                            }
                        }
                    }
                },
                "Superconduction temperature": "0.2 K",
                "Thermal conductivity": "27 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.86,
                "Velocity of sound": "3155 m s<sup>-1</sup>",
                "Vickers hardness": "1960 MN m<sup>-2</sup>",
                "X": 1.38,
                "Youngs modulus": "208 GPa",
                "Metallic radius": 1.542,
                "iupac_ordering": 29,
                "IUPAC ordering": 29
            },
            "V": {
                "Atomic mass": 50.9415,
                "Atomic no": 23,
                "Atomic orbitals": {
                    "1s": -195.224014,
                    "2p": -18.435189,
                    "2s": -21.815346,
                    "3d": -0.204634,
                    "3p": -1.610516,
                    "3s": -2.526904,
                    "4s": -0.175968
                },
                "Atomic radius": 1.35,
                "Atomic radius calculated": 1.71,
                "Boiling point": "3680 K",
                "Brinell hardness": "628 MN m<sup>-2</sup>",
                "Bulk modulus": "160 GPa",
                "Coefficient of linear thermal expansion": "8.4 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    5
                ],
                "Critical temperature": "no data K",
                "Density of solid": "6110 kg m<sup>-3</sup>",
                "Electrical resistivity": "20 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ar].3d<sup>3</sup>.4s<sup>2</sup>",
                "ICSD oxidation states": [
                    2,
                    3,
                    4,
                    5
                ],
                "Ionic radii": {
                    "2": 0.93,
                    "3": 0.78,
                    "4": 0.72,
                    "5": 0.68
                },
                "Liquid range": "1497 K",
                "Melting point": "2183 K",
                "Mendeleev no": 54,
                "Mineral hardness": "7.0",
                "Molar volume": "8.32 cm<sup>3</sup>",
                "Name": "Vanadium",
                "Oxidation states": [
                    -1,
                    1,
                    2,
                    3,
                    4,
                    5
                ],
                "Poissons ratio": "0.37",
                "Reflectivity": "61 %",
                "Refractive index": "no data",
                "Rigidity modulus": "47 GPa",
                "Shannon radii": {
                    "2": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.93,
                                "ionic_radius": 0.79
                            }
                        }
                    },
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.78,
                                "ionic_radius": 0.64
                            }
                        }
                    },
                    "4": {
                        "V": {
                            "": {
                                "crystal_radius": 0.67,
                                "ionic_radius": 0.53
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.72,
                                "ionic_radius": 0.58
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 0.86,
                                "ionic_radius": 0.72
                            }
                        }
                    },
                    "5": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.495,
                                "ionic_radius": 0.355
                            }
                        },
                        "V": {
                            "": {
                                "crystal_radius": 0.6,
                                "ionic_radius": 0.46
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.68,
                                "ionic_radius": 0.54
                            }
                        }
                    }
                },
                "Superconduction temperature": "5.40 K",
                "Thermal conductivity": "31 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "4560 m s<sup>-1</sup>",
                "Vickers hardness": "628 MN m<sup>-2</sup>",
                "X": 1.63,
                "Youngs modulus": "128 GPa",
                "NMR Quadrupole Moment": {
                    "V-50": 210.4,
                    "V-51": -52.1
                },
                "Metallic radius": 1.347,
                "iupac_ordering": 55,
                "IUPAC ordering": 55
            },
            "W": {
                "Atomic mass": 183.84,
                "Atomic no": 74,
                "Atomic orbitals": {
                    "1s": -2341.042887,
                    "2p": -369.013973,
                    "2s": -384.856157,
                    "3d": -66.724787,
                    "3p": -80.502102,
                    "3s": -87.867792,
                    "4d": -8.879693,
                    "4f": -1.550835,
                    "4p": -14.495102,
                    "4s": -17.570797,
                    "5d": -0.220603,
                    "5p": -1.504457,
                    "5s": -2.396018,
                    "6s": -0.181413
                },
                "Atomic radius": 1.35,
                "Atomic radius calculated": 1.93,
                "Boiling point": "5828 K",
                "Brinell hardness": "2570 MN m<sup>-2</sup>",
                "Bulk modulus": "310 GPa",
                "Coefficient of linear thermal expansion": "4.5 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    4,
                    6
                ],
                "Critical temperature": "no data K",
                "Density of solid": "19250 kg m<sup>-3</sup>",
                "Electrical resistivity": "5.4 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>14</sup>.5d<sup>4</sup>.6s<sup>2</sup>",
                "ICSD oxidation states": [
                    2,
                    3,
                    4,
                    5,
                    6
                ],
                "Ionic radii": {
                    "4": 0.8,
                    "5": 0.76,
                    "6": 0.74
                },
                "Liquid range": "2133 K",
                "Melting point": "3695 K",
                "Mendeleev no": 55,
                "Mineral hardness": "7.5",
                "Molar volume": "9.47 cm<sup>3</sup>",
                "Name": "Tungsten",
                "Oxidation states": [
                    -2,
                    -1,
                    1,
                    2,
                    3,
                    4,
                    5,
                    6
                ],
                "Poissons ratio": "0.28",
                "Reflectivity": "62 %",
                "Refractive index": "no data",
                "Rigidity modulus": "161 GPa",
                "Shannon radii": {
                    "4": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.8,
                                "ionic_radius": 0.66
                            }
                        }
                    },
                    "5": {
                        "VI": {
                            "": {
                                "crystal_radius": 0.76,
                                "ionic_radius": 0.62
                            }
                        }
                    },
                    "6": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.56,
                                "ionic_radius": 0.42
                            }
                        },
                        "V": {
                            "": {
                                "crystal_radius": 0.65,
                                "ionic_radius": 0.51
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.74,
                                "ionic_radius": 0.6
                            }
                        }
                    }
                },
                "Superconduction temperature": "0.015 K",
                "Thermal conductivity": "170 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "5174 m s<sup>-1</sup>",
                "Vickers hardness": "3430 MN m<sup>-2</sup>",
                "X": 2.36,
                "Youngs modulus": "411 GPa",
                "Metallic radius": 1.41,
                "iupac_ordering": 56,
                "IUPAC ordering": 56
            },
            "Xe": {
                "Atomic mass": 131.293,
                "Atomic no": 54,
                "Atomic orbitals": {
                    "1s": -1208.688993,
                    "2p": -172.599583,
                    "2s": -183.327495,
                    "3d": -24.37823,
                    "3p": -32.867042,
                    "3s": -37.415454,
                    "4d": -2.286666,
                    "4p": -5.063802,
                    "4s": -6.67834,
                    "5p": -0.309835,
                    "5s": -0.672086
                },
                "Atomic radius": "no data",
                "Atomic radius calculated": 1.08,
                "Boiling point": "165.1 K",
                "Brinell hardness": "no data MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "no data x10<sup>-6</sup>K<sup>-1</sup>",
                "Critical temperature": "289.7 K",
                "Density of solid": "no data kg m<sup>-3</sup>",
                "Electrical resistivity": "no data 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Kr].4d<sup>10</sup>.5s<sup>2</sup>.5p<sup>6</sup>",
                "Ionic radii": {
                    "8": 0.62
                },
                "Liquid range": "3.7 K",
                "Max oxidation state": 8.0,
                "Melting point": "161.4 K",
                "Mendeleev no": 5,
                "Min oxidation state": 2.0,
                "Mineral hardness": "no data",
                "Molar volume": "35.92 cm<sup>3</sup>",
                "Name": "Xenon",
                "Poissons ratio": "no data",
                "Reflectivity": "no data %",
                "Refractive index": "1.000702",
                "Rigidity modulus": "no data GPa",
                "Shannon radii": {
                    "8": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.54,
                                "ionic_radius": 0.4
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.62,
                                "ionic_radius": 0.48
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "0.00565 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 2.16,
                "Velocity of sound": "1090 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 2.6,
                "Youngs modulus": "no data GPa",
                "Metallic radius": "no data",
                "iupac_ordering": 1,
                "IUPAC ordering": 1
            },
            "Y": {
                "Atomic mass": 88.90585,
                "Atomic no": 39,
                "Atomic orbitals": {
                    "1s": -605.631981,
                    "2p": -74.803201,
                    "2s": -81.789102,
                    "3d": -5.671499,
                    "3p": -10.399926,
                    "3s": -12.992217,
                    "4d": -0.108691,
                    "4p": -1.02449,
                    "4s": -1.697124,
                    "5s": -0.150727
                },
                "Atomic radius": 1.8,
                "Atomic radius calculated": 2.12,
                "Boiling point": "3609 K",
                "Brinell hardness": "589 MN m<sup>-2</sup>",
                "Bulk modulus": "41 GPa",
                "Coefficient of linear thermal expansion": "10.6 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "4472 kg m<sup>-3</sup>",
                "Electrical resistivity": "about 60 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Kr].4d<sup>1</sup>.5s<sup>2</sup>",
                "ICSD oxidation states": [
                    3
                ],
                "Ionic radii": {
                    "3": 1.04
                },
                "Liquid range": "1810 K",
                "Melting point": "1799 K",
                "Mendeleev no": 25,
                "Mineral hardness": "no data",
                "Molar volume": "19.88 cm<sup>3</sup>",
                "Name": "Yttrium",
                "Oxidation states": [
                    1,
                    2,
                    3
                ],
                "Poissons ratio": "0.24",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "26 GPa",
                "Shannon radii": {
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.04,
                                "ionic_radius": 0.9
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.1,
                                "ionic_radius": 0.96
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.159,
                                "ionic_radius": 1.019
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.215,
                                "ionic_radius": 1.075
                            }
                        }
                    }
                },
                "Superconduction temperature": "1.3  (under pressure)K",
                "Thermal conductivity": "17 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "3300 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.22,
                "Youngs modulus": "64 GPa",
                "Metallic radius": 1.8,
                "iupac_ordering": 48,
                "IUPAC ordering": 48
            },
            "Yb": {
                "Atomic mass": 173.04,
                "Atomic no": 70,
                "Atomic orbitals": {
                    "1s": -2084.069389,
                    "2p": -323.178219,
                    "2s": -337.978976,
                    "3d": -56.026315,
                    "3p": -68.698655,
                    "3s": -75.47663,
                    "4d": -6.574963,
                    "4f": -0.286408,
                    "4p": -11.558246,
                    "4s": -14.312076,
                    "5p": -0.966137,
                    "5s": -1.683886,
                    "6s": -0.136989
                },
                "Atomic radius": 1.75,
                "Atomic radius calculated": 2.22,
                "Boiling point": "1469 K",
                "Brinell hardness": "343 MN m<sup>-2</sup>",
                "Bulk modulus": "31 GPa",
                "Coefficient of linear thermal expansion": "26.3 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    3
                ],
                "Critical temperature": "no data K",
                "Density of solid": "6570 kg m<sup>-3</sup>",
                "Electrical resistivity": "25.0 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Xe].4f<sup>14</sup>.6s<sup>2</sup>",
                "ICSD oxidation states": [
                    2,
                    3
                ],
                "Ionic radii": {
                    "2": 1.16,
                    "3": 1.008
                },
                "Liquid range": "372 K",
                "Melting point": "1097 K",
                "Mendeleev no": 17,
                "Mineral hardness": "no data",
                "Molar volume": "24.84 cm<sup>3</sup>",
                "Name": "Ytterbium",
                "Oxidation states": [
                    2,
                    3
                ],
                "Poissons ratio": "0.21",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "9.9 GPa",
                "Shannon radii": {
                    "2": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.16,
                                "ionic_radius": 1.02
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.22,
                                "ionic_radius": 1.08
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.28,
                                "ionic_radius": 1.14
                            }
                        }
                    },
                    "3": {
                        "VI": {
                            "": {
                                "crystal_radius": 1.008,
                                "ionic_radius": 0.868
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 1.065,
                                "ionic_radius": 0.925
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.125,
                                "ionic_radius": 0.985
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.182,
                                "ionic_radius": 1.042
                            }
                        }
                    }
                },
                "Superconduction temperature": "no data K",
                "Thermal conductivity": "39 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "1590 m s<sup>-1</sup>",
                "Vickers hardness": "206 MN m<sup>-2</sup>",
                "X": 1.1,
                "Youngs modulus": "24 GPa",
                "Metallic radius": 1.94,
                "iupac_ordering": 34,
                "IUPAC ordering": 34
            },
            "Zn": {
                "Atomic mass": 65.409,
                "Atomic no": 30,
                "Atomic orbitals": {
                    "1s": -344.969756,
                    "2p": -36.648765,
                    "2s": -41.531323,
                    "3d": -0.398944,
                    "3p": -3.022363,
                    "3s": -4.573041,
                    "4s": -0.222725
                },
                "Atomic radius": 1.35,
                "Atomic radius calculated": 1.42,
                "Boiling point": "1180 K",
                "Brinell hardness": "412 MN m<sup>-2</sup>",
                "Bulk modulus": "70 GPa",
                "Coefficient of linear thermal expansion": "30.2 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    2
                ],
                "Critical temperature": "no data K",
                "Density of solid": "7140 kg m<sup>-3</sup>",
                "Electrical resistivity": "6.0 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Ar].3d<sup>10</sup>.4s<sup>2</sup>",
                "ICSD oxidation states": [
                    2
                ],
                "Ionic radii": {
                    "2": 0.88
                },
                "Liquid range": "487.32 K",
                "Melting point": "692.68 K",
                "Mendeleev no": 76,
                "Mineral hardness": "2.5",
                "Molar volume": "9.16 cm<sup>3</sup>",
                "Name": "Zinc",
                "Oxidation states": [
                    1,
                    2
                ],
                "Poissons ratio": "0.25",
                "Reflectivity": "80 %",
                "Refractive index": "1.002050",
                "Rigidity modulus": "43 GPa",
                "Shannon radii": {
                    "2": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.74,
                                "ionic_radius": 0.6
                            }
                        },
                        "V": {
                            "": {
                                "crystal_radius": 0.82,
                                "ionic_radius": 0.68
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.88,
                                "ionic_radius": 0.74
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 1.04,
                                "ionic_radius": 0.9
                            }
                        }
                    }
                },
                "Superconduction temperature": "0.85 K",
                "Thermal conductivity": "120 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": 1.39,
                "Velocity of sound": "3700 m s<sup>-1</sup>",
                "Vickers hardness": "no data MN m<sup>-2</sup>",
                "X": 1.65,
                "Youngs modulus": "108 GPa",
                "NMR Quadrupole Moment": {
                    "Zn-67": 150.15
                },
                "Metallic radius": 1.34,
                "iupac_ordering": 76,
                "IUPAC ordering": 76
            },
            "Zr": {
                "Atomic mass": 91.224,
                "Atomic no": 40,
                "Atomic orbitals": {
                    "1s": -639.292236,
                    "2p": -80.010043,
                    "2s": -87.237062,
                    "3d": -6.544643,
                    "3p": -11.514415,
                    "3s": -14.230432,
                    "4d": -0.150673,
                    "4p": -1.186597,
                    "4s": -1.918971,
                    "5s": -0.162391
                },
                "Atomic radius": 1.55,
                "Atomic radius calculated": 2.06,
                "Boiling point": "4682 K",
                "Brinell hardness": "650 MN m<sup>-2</sup>",
                "Bulk modulus": "no data GPa",
                "Coefficient of linear thermal expansion": "5.7 x10<sup>-6</sup>K<sup>-1</sup>",
                "Common oxidation states": [
                    4
                ],
                "Critical temperature": "no data K",
                "Density of solid": "6511 kg m<sup>-3</sup>",
                "Electrical resistivity": "43.3 10<sup>-8</sup> &Omega; m",
                "Electronic structure": "[Kr].4d<sup>2</sup>.5s<sup>2</sup>",
                "ICSD oxidation states": [
                    2,
                    3,
                    4
                ],
                "Ionic radii": {
                    "4": 0.86
                },
                "Liquid range": "2554 K",
                "Melting point": "2128 K",
                "Mendeleev no": 49,
                "Mineral hardness": "5.0",
                "Molar volume": "14.02 cm<sup>3</sup>",
                "Name": "Zirconium",
                "Oxidation states": [
                    1,
                    2,
                    3,
                    4
                ],
                "Poissons ratio": "0.34",
                "Reflectivity": "no data %",
                "Refractive index": "no data",
                "Rigidity modulus": "33 GPa",
                "Shannon radii": {
                    "4": {
                        "IV": {
                            "": {
                                "crystal_radius": 0.73,
                                "ionic_radius": 0.59
                            }
                        },
                        "V": {
                            "": {
                                "crystal_radius": 0.8,
                                "ionic_radius": 0.66
                            }
                        },
                        "VI": {
                            "": {
                                "crystal_radius": 0.86,
                                "ionic_radius": 0.72
                            }
                        },
                        "VII": {
                            "": {
                                "crystal_radius": 0.92,
                                "ionic_radius": 0.78
                            }
                        },
                        "VIII": {
                            "": {
                                "crystal_radius": 0.98,
                                "ionic_radius": 0.84
                            }
                        },
                        "IX": {
                            "": {
                                "crystal_radius": 1.03,
                                "ionic_radius": 0.89
                            }
                        }
                    }
                },
                "Superconduction temperature": "0.61 K",
                "Thermal conductivity": "23 W m<sup>-1</sup> K<sup>-1</sup>",
                "Van der waals radius": "no data",
                "Velocity of sound": "3800 m s<sup>-1</sup>",
                "Vickers hardness": "903 MN m<sup>-2</sup>",
                "X": 1.33,
                "Youngs modulus": "68 GPa",
                "Metallic radius": 1.602,
                "iupac_ordering": 51,
                "IUPAC ordering": 51
            },
            "Rf": {
                "Atomic mass": 267,
                "Atomic no": 104,
                "Name": "Rutherfordium"
            },
            "Db": {
                "Atomic mass": 268,
                "Atomic no": 105,
                "Name": "Dubnium"
            },
            "Sg": {
                "Atomic mass": 269,
                "Atomic no": 106,
                "Name": "Seaborgium"
            },
            "Bh": {
                "Atomic mass": 270,
                "Atomic no": 107,
                "Name": "Bohrium"
            },
            "Hs": {
                "Atomic mass": 270,
                "Atomic no": 108,
                "Name": "Hassium"
            },
            "Mt": {
                "Atomic mass": 278,
                "Atomic no": 109,
                "Name": "Meitnerium"
            },
            "Ds": {
                "Atomic mass": 281,
                "Atomic no": 110,
                "Name": "Darmstadtium"
            },
            "Rg": {
                "Atomic mass": 282,
                "Atomic no": 111,
                "Name": "Roentgenium"
            },
            "Cn": {
                "Atomic mass": 285,
                "Atomic no": 112,
                "Name": "Copernicium"
            },
            "Nh": {
                "Atomic mass": 286,
                "Atomic no": 113,
                "Name": "Nihonium"
            },
            "Fl": {
                "Atomic mass": 289,
                "Atomic no": 114,
                "Name": "Flerovium"
            },
            "Mc": {
                "Atomic mass": 290,
                "Atomic no": 115,
                "Name": "Moscovium"
            },
            "Lv": {
                "Atomic mass": 293,
                "Atomic no": 116,
                "Name": "Livermorium"
            },
            "Ts": {
                "Atomic mass": 294,
                "Atomic no": 117,
                "Name": "Tennessine"
            },
            "Og": {
                "Atomic mass": 2949,
                "Atomic no": 118,
                "Name": "Oganesson"
            }
        }

        self.enmax = {'Ac': 172.351,
                      'Ag': 249.844,
                      'Al': 240.3,
                      'Am': 255.875,
                      'Ar': 266.408,
                      'As': 208.702,
                      'At': 161.43,
                      'Au': 229.943,
                      'B': 318.614,
                      'Ba': 187.181,
                      'Be': 247.543,
                      'Bi': 105.037,
                      'Br': 216.285,
                      'C': 400.0,
                      'Ca': 266.622,
                      'Cd': 274.336,
                      'Ce': 273.042,
                      'Cf': 414.614,
                      'Cl': 262.472,
                      'Cm': 257.953,
                      'Co': 267.968,
                      'Cr': 227.08,
                      'Cs': 220.318,
                      'Cu': 295.446,
                      'Dy': 255.467,
                      'Er': 298.116,
                      'Eu': 249.668,
                      'F': 400.0,
                      'Fe': 267.882,
                      'Fr': 214.54,
                      'Ga': 134.678,
                      'Gd': 256.472,
                      'Ge': 173.807,
                      'H': 250.0,
                      'He': 478.896,
                      'Hf': 220.334,
                      'Hg': 233.204,
                      'Ho': 257.168,
                      'I': 175.647,
                      'In': 95.934,
                      'Ir': 210.864,
                      'Kr': 185.331,
                      'K': 259.264,
                      'La': 219.292,
                      'Li': 140.0,
                      'Lu': 255.695,
                      'Mg': 200.0,
                      'Mn': 269.864,
                      'Mo': 224.584,
                      'N': 400.0,
                      'Na': 101.968,
                      'Nb': 293.235,
                      'Nd': 253.189,
                      'Ne': 343.606,
                      'Ni': 269.532,
                      'Np': 254.26,
                      'O': 400.0,
                      'Os': 228.022,
                      'P': 255.04,
                      'Pa': 252.193,
                      'Pb': 97.973,
                      'Pd': 250.925,
                      'Pm': 258.627,
                      'Po': 159.707,
                      'Pr': 272.941,
                      'Pt': 230.283,
                      'Pu': 254.353,
                      'Ra': 237.367,
                      'Rb': 220.112,
                      'Re': 226.216,
                      'Rh': 228.996,
                      'Rn': 152.121,
                      'Ru': 213.271,
                      'S': 258.689,
                      'Sb': 172.069,
                      'Sc': 154.763,
                      'Se': 211.555,
                      'Si': 245.345,
                      'Sm': 257.515,
                      'Sn': 103.236,
                      'Sr': 229.353,
                      'Ta': 223.667,
                      'Tb': 264.824,
                      'Tc': 228.694,
                      'Te': 174.982,
                      'Th': 247.306,
                      'Ti': 178.33,
                      'Tl': 90.14,
                      'Tm': 257.42,
                      'U': 252.502,
                      'V': 192.543,
                      'W': 223.057,
                      'Xe': 153.118,
                      'Yb': 253.028,
                      'Y': 202.626,
                      'Zn': 276.723,
                      'Zr': 229.898}

        self.atom_energy = {
            # 2015npj_OQMD assessing the accuracy of DFT formation energies
            'H': -3.434,
            'He': -0.004,
            'Li': -1.731,
            'Be': -3.653,
            "B": -6.656,
            "C": -9.044,
            "N": -8.195,
            "O": -4.523,
            "F": -1.443,
            'Ne': -0.029,
            'Na': -1.196,
            'Mg': -1.417,
            'Al': -3.660,
            'Si': -5.386,
            'P': -5.175,
            'S': -3.868,
            'Cl': -1.479,
            'Ar': -0.006,
            'K': -0.987,
            'Ca': -1.780,
            'Sc': -6.344,
            'Ti': -7.702,
            'V': -8.898,
            'Cr': -9.463,
            'Mn': -8.898,
            'Fe': -8.499,
            'Co': -7.078,
            'Ni': -5.587,
            'Cu': -3.710,
            'Zn': -1.157,
            'Ga': -2.902,
            'Ge': -4.522,
            'As': -4.593,
            'Se': -3.374,
            'Br': -1.333,
            'Kr': -0.004,
            'Rb': -0.881,
            'Sr': -1.549,
            'Y': -6.449,
            'Zr': -8.438,
            'Nb': -10.017,
            'Mo': -10.921,
            'Tc': -10.457,
            'Ru': -9.210,
            'Rh': -7.319,
            'Pd': -5.197,
            'Ag': -2.907,
            'Cd': -0.861,
            'In': -2.609,
            'Sn': -3.938,
            'Sb': -4.155,
            'Te': -3.027,
            'I': -1.365,
            'Xe': -0.640,
            'Cs': -0.744,
            'Ba': -1.479,
            'La': -4.959,
            'Ce': -4.564,
            'Pr': -4.627,
            'Nd': -4.697,
            'Pm': -4.716,
            'Sm': -4.606,
            'Eu': -1.708,
            'Gd': -4.718,
            'Tb': -4.731,
            'Dy': -4.660,
            'Ho': -4.573,
            'Er': -4.582,
            'Tm': -4.451,
            'Yb': -1.125,
            'Lu': -4.549,
            'Hf': -9.902,
            'Ta': -11.941,
            'W': -13.130,
            'Re': -12.378,
            'Os': -11.374,
            'Ir': -8.953,
            'Pt': -6.162,
            'Au': -3.283,
            'Hg': -0.374,
            'Tl': -2.480,
            'Pb': -3.951,
            'Bi': -4.199,
            'Ac': -4.106,
            'Th': -7.237,
            'Pa': -9.497,
            'U': -11.032,
            'Np': -12.797,
            'Pu': -13.950
        }

