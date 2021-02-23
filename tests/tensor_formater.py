#!/usr/bin/python
import itertools
import sys

def print_with_braces(shape, data, idx, depth = 0):
    """
    This fxn handles the recursion necessary to print a tensor out with the
    proper number of braces.
    """
    rank = len(shape)
    if depth == rank:
        return ''
    rv = '{'
    for i in range(shape[depth]):
        idx[depth] = i
        if depth == rank - 1:
            rv += data[tuple(idx)]
        else:
            rv += print_with_braces(shape, data,  idx, depth + 1)
        if i < shape[depth] - 1:
            rv += ', '
            if depth != rank - 1:
                rv += '\n'
    rv += '}'
    return rv        


def ta_output_to_input(file_2_read_from):
    """ This function assumes that file_2_read_from contains the output of
        printing a TiledArray tensor.
    """    
    with open(file_2_read_from, 'r') as f:
        data = {}
        shape = None
        for line in f.readlines():
            if not '{' in line:
                continue
            tile_start = line.split('[')[2].split(']')[0].split(',')
            tile_end   = line.split('[')[3].split(']')[0].split(',')
            str_data   = line.split('{')[1].split('}')[0].split()
            
            if shape:
                shape = [max(shape[i], int(te)) for i, te in enumerate(tile_end)]
            else:
                shape = [int(te) for te in tile_end]

            ranges = []
            for mode_start, mode_end in zip(tile_start, tile_end):
                ranges.append(range(int(mode_start), int(mode_end)))

            for i, idx in enumerate(itertools.product(*ranges)):
                data[idx] = str_data[i]

        zero_idx = [0 for _ in range(len(shape))]
        rv = print_with_braces(shape, data, zero_idx)
        print(rv)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        msg = 'Syntax: tensor_formater.py file_with_tensor'
        raise Exception(msg)
    ta_output_to_input(sys.argv[1])
