#!/usr/local/bin/python

import argparse
from collections import OrderedDict
import numpy as np
import math


def vectorlength(vec):
    len_vec = math.sqrt(math.pow(float(vec[0]), 2.0) + math.pow(float(vec[1]), 2.0) + math.pow(float(vec[2]), 2.0))
    return len_vec


def dot(vec1, vec2):
    dot_prod = 0.0

    for i in range(3):
        dot_prod += (float(vec1[i]) * float(vec2[i]))

    return dot_prod


def angle(vec1, vec2):
    angle = math.degrees(math.acos(round(dot(vec1, vec2) /
                                         (vectorlength(vec1) * vectorlength(vec2)), 10)))
    return angle


def cellmatrix(unitvec):
    alpha = angle(unitvec[1], unitvec[2])
    beta = angle(unitvec[2], unitvec[0])
    gamma = angle(unitvec[0], unitvec[1])

    ang = [alpha, beta, gamma]
    angles = np.array(ang)

    a = vectorlength(unitvec[0])
    b = vectorlength(unitvec[1])
    c = vectorlength(unitvec[2])

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

    return fileindex, unitvec, atominfo, cart, seldyn, coord, dyn, vel, atomlist


def distancecalc(cellinfo, output, species, tolerance):

    return


def anglecalc(cellinfo, output, species, tolerance):
    return


def calculatedistance(args):
    p = parser(args.input)
    d = distancecalc(p, args.output, args.species, args.tolerance)
    return


def calculateangle(args):
    p = parser(args.input)
    d = anglecalc(p, args.output, args.species, args.tolerance)
    return


def main():
    # Main parser
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)

    # Subparsers
    subparsers = parser.add_subparsers(title="Functions")

    parser_dist = subparsers.add_parser("dist", formatter_class=argparse.RawTextHelpFormatter, description=descpoly)
    parser_dist.add_argument("-i", dest="input", type=str, default="POSCAR")
    parser_dist.add_argument("-o", dest="output", type=str, default=None)
    parser_dist.add_argument("-s", dest="species", type=str, required=True, nargs="*")
    parser_dist.add_argument("-t", dest="tol", type=float, default=0.01)
    parser_dist.set_defaults(func=calculatedistance)

    parser_angle = subparsers.add_parser("angle", formatter_class=argparse.RawTextHelpFormatter, description=descpoly)
    parser_angle.add_argument("-i", dest="input", type=str, default="POSCAR")
    parser_angle.add_argument("-o", dest="output", type=str, default=None)
    parser_angle.add_argument("-s", dest="species", type=str, required=True, nargs="*")
    parser_angle.add_argument("-t", dest="tol", type=float, default=0.1)
    parser_angle.set_defaults(func=calculateangle)


