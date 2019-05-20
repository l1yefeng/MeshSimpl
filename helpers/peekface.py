#!/usr/bin/env python3
#
# Print lines in OBJ file with face index (0-based)

from sys import argv

if __name__ == '__main__':
    if len(argv) < 3:
        print('Usage: {} OBJFILE FID [FID...]'.format(argv[0]))
        exit(1)

    filepath = argv[1]
    targets = sorted(int(f) for f in argv[2:])

    pos = 0
    idx = 0
    with open(filepath, 'r') as f:
        for line in f:
            if not line.startswith('f '):
                continue
            if idx == targets[pos]:
                print('f#{} -> {}'.format(idx, line.strip()))
                pos += 1
                if pos == len(targets):
                    break
            idx += 1

    for f in targets[pos:]:
        print('f#{} -> not in range'.format(f))
