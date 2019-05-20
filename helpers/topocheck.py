#!/usr/bin/env python3
#
# Find non-manifold edges in provided mesh

from sys import argv
from typing import List, Dict, Tuple


def parse_obj(filepath: str) -> List[List[int]]:
    indices = []
    with open(filepath, 'r') as f:
        for line in f:
            if not line.startswith('f '):
                continue
            toks = line.split()
            assert (len(toks) == 4)
            toks = toks[1:]
            corners = [int(x.partition('/')[0]) for x in toks]
            indices.append(corners)
    return indices


if __name__ == '__main__':
    if len(argv) != 2:
        print('Usage: {} OBJFILE'.format(argv[0]))
        exit(1)

    indices = parse_obj(argv[1])

    edges: Dict[Tuple[int, int], List[int]] = dict()

    for f, face in enumerate(indices):
        for i, j in [(1, 2), (2, 0), (0, 1)]:
            v0, v1 = face[i], face[j]
            assert (v0 != v1)
            if v0 > v1:
                v0, v1 = v1, v0
            key = (v0, v1)
            efaces = edges.get(key)
            if efaces is None:
                edges[key] = [f]
            else:
                efaces.append(f)

    nonmani = [faces for _, faces in edges.items() if len(faces) > 2]

    cnt = len(nonmani)
    if cnt == 0:
        print('OK')
    else:
        if cnt == 1:
            print('1 non-manifold edge was found')
        else:
            print('{} non-manifold edges were found'.format(cnt))
        for efaces in nonmani:
            for i, f in enumerate(efaces):
                line = '->' if i == 0 else '  '
                line += ' f{}:'.format(f)
                for v in indices[f]:
                    line += '\t' + str(v + 1)
                print(line)
