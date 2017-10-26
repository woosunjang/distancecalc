#!/usr/local/bin/python

import argparse
import numpy as np
import math
import itertools
from collections import OrderedDict


def vectorlength(vec):
    len_vec = math.sqrt(math.pow(float(vec[0]), 2.0) + math.pow(float(vec[1]), 2.0) + math.pow(float(vec[2]), 2.0))
    return len_vec


def vectordistance(vec1, vec2):
    dist_vec = math.sqrt(math.pow(float(vec1[0]) - float(vec2[0]), 2.0) +
                         math.pow(float(vec1[1]) - float(vec2[1]), 2.0) +
                         math.pow(float(vec1[2]) - float(vec2[2]), 2.0))
    return dist_vec


def dot(vec1, vec2):
    dot_prod = 0.0

    for i in range(3):
        dot_prod += (float(vec1[i]) * float(vec2[i]))

    return dot_prod


def angle(vec1, vec2):
    ang = math.degrees(math.acos(round(dot(vec1, vec2) / (vectorlength(vec1) * vectorlength(vec2)), 10)))
    return ang


def cellmatrix(unitvec):
    alpha = angle(unitvec[1], unitvec[2])
    beta = angle(unitvec[2], unitvec[0])
    gamma = angle(unitvec[0], unitvec[1])

    ang = [alpha, beta, gamma]
    angles = np.array(ang)

    mat = unitvec
    matrix = np.array(mat)
    return matrix, angles


def parser(infile):
    print("Reading " + str(infile) + "...")

    with open(infile, 'r') as f:
        fileindex = f.readline().strip(' \n')
        scalfac = float(f.readline().strip(' \n'))

        unitvec = []
        for i in range(3):
            unitvec.append(f.readline().split())
        unitvec = np.array(unitvec, dtype='d') * scalfac

        atom_species = f.readline().split()
        atom_counts = f.readline().split()
        atominfo = OrderedDict()
        for i in range(len(atom_species)):
            atominfo[atom_species[i]] = int(atom_counts[i])
        atominfo = atominfo
        count = 0
        for x in atominfo:
            count += atominfo[x]

        determ_line = f.readline()
        if str(determ_line[0]) == "S":
            seldyn = True
            coordtype = f.readline().strip(' \n')
        else:
            seldyn = False
            coordtype = determ_line.strip(' \n')

        cart = None
        if coordtype[0] == 'D':
            cart = False
        elif coordtype[0] == 'd':
            cart = False
        elif coordtype[0] == 'C':
            cart = True
        elif coordtype[0] == 'c':
            cart = True

        tmp = f.read()

        coord = []
        dyn = []
        vel = []

        if not seldyn:
            tmp = np.array(tmp.split(), dtype='d')
            tmp = np.reshape(tmp, (int(len(tmp) / 3), 3))
            dyn = np.array([], dtype='str')
            for i in range(count):
                coord.append(tmp[i])
                if len(tmp) == count * 2:
                    vel.append(tmp[i + count])
            coord = np.array(coord)
            vel = np.array(vel)

        elif seldyn:
            tmp = np.array(tmp.split())
            tmpc = tmp[:(count * 6)]
            if len(tmp) == count * 9:
                vel = np.reshape(tmp[(count * 6):], (count, 3))
            vel = np.array(vel, dtype='d')
            for i in range(count):
                for j in range(3):
                    coord.append(tmpc[i * 6 + j])
                    dyn.append(tmpc[i * 6 + j + 3])
            coord = np.reshape(np.array(coord, dtype='d'), (count, 3))
            dyn = np.reshape(np.array(dyn, dtype='str'), (count, 3))

        atomlist = []
        for i in atominfo:
            for j in range(atominfo[i]):
                atomlist.append(i)
        atomlist = np.vstack(atomlist)

        dic = {"index": fileindex,
               "unitvec": unitvec,
               "atominfo": atominfo,
               "cart": cart,
               "seldyn": seldyn,
               "coord": coord,
               "dyn": dyn,
               "vel": vel,
               "atomlist": atomlist
               }

    return dic


def distancecalc(cellinfo, output, species, tolerance):
    dist = []
    comb = []
    ind = []

    if cellinfo['cart'] is False:
        pos = np.dot(cellinfo['coord'], cellinfo['unitvec'])
    else:
        pos = cellinfo['coord']

    if len(species) == 1:
        for x in itertools.combinations_with_replacement(species, 2):
            comb.append(x)

    else:
        for x in itertools.combinations(species, 2):
            comb.append(x)

    for x in comb:
        for y in itertools.combinations(enumerate(cellinfo['atomlist']), 2):
            for z in itertools.combinations([y[0][1][0], y[1][1][0]], 2):
                if x == z:
                    ind.append([y[0][0], y[1][0]])

    for x in ind:
        dist.append(vectordistance(pos[x[0]], pos[x[1]]))

    return np.array(sorted(dist, key=float))


def anglecalc(cellinfo, output, species, tolerance):
    ang = []
    comb = []
    ind = []

    if cellinfo['cart'] is False:
        pos = np.dot(cellinfo['coord'], cellinfo['unitvec'])
    else:
        pos = cellinfo['coord']

    if len(species) == 1:
        for x in itertools.combinations_with_replacement(species, 3):
            comb.append(x)

    else:
        for x in itertools.combinations(species, 3):
            comb.append(x)

    for x in comb:
        for y in itertools.combinations(enumerate(cellinfo['atomlist']), 3):
            for z in itertools.combinations([y[0][1][0], y[1][1][0], y[2][1][0]], 3):
                if x == z:
                    ind.append([y[0][0], y[1][0], y[2][0]])

    for x in ind:
        vec1 = pos[x[0]] - pos[x[2]]
        vec2 = pos[x[1]] - pos[x[2]]
        ang.append(angle(vec1, vec2))

    return np.array(sorted(ang, key=float))


def calculatedistance(args):
    p = parser(args.input)
    d = distancecalc(p, args.output, args.species, args.tol)
    return


def calculateangle(args):
    p = parser(args.input)
    d = anglecalc(p, args.output, args.species, args.tol)
    return


def main():
    description = ""
    descdist = ""
    descangle = ""

    # Main parser
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)

    # Subparsers
    subparsers = parser.add_subparsers(title="Functions")

    parser_dist = subparsers.add_parser("dist", formatter_class=argparse.RawTextHelpFormatter, description=descdist)
    parser_dist.add_argument("-i", dest="input", type=str, default="POSCAR")
    parser_dist.add_argument("-o", dest="output", type=str, default=None)
    parser_dist.add_argument("-s", dest="species", type=str, required=True, nargs="*")
    parser_dist.add_argument("-t", dest="tol", type=float, default=0.01)
    parser_dist.set_defaults(func=calculatedistance)

    parser_angle = subparsers.add_parser("angle", formatter_class=argparse.RawTextHelpFormatter, description=descangle)
    parser_angle.add_argument("-i", dest="input", type=str, default="POSCAR")
    parser_angle.add_argument("-o", dest="output", type=str, default=None)
    parser_angle.add_argument("-s", dest="species", type=str, required=True, nargs="*")
    parser_angle.add_argument("-t", dest="tol", type=float, default=0.1)
    parser_angle.set_defaults(func=calculateangle)

    args = parser.parse_args()

    try:
        getattr(args, "func")
    except AttributeError:
        parser.print_help()
        raise AttributeError("Error!")
    args.func(args)

if __name__ == "__main__":
    main()