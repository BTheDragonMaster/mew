def write_bpp(well_to_bpp, out_file):
    with open(out_file, 'w') as out:
        for well, bpps in well_to_bpp.items():
            out.write(f'{well}')
            for bpp in bpps:
                out.write(f',{bpp:.10f}')
            out.write('\n')


def write_sequences(well_to_sequence, out_file):
    with open(out_file, 'w') as out:
        for well, sequence in well_to_sequence.items():
            out.write(f'{well},{sequence}\n')


def write_flow(well_to_flow, out_file):
    with open(out_file, 'w') as out:
        for well, flow in well_to_flow.items():
            out.write(f'{well},{flow}\n')