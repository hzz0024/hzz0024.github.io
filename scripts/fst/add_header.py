import sys
with open(sys.argv[1], 'r') as fin:
    lines = fin.readlines() 
    header = lines[0]
    _, t1, t2 = header.strip('\n').split()
    lines[0] = lines[0].replace(t1, 'window_start').replace(t2, 'window_end angsd_Fst')
with open(sys.argv[2], 'w') as fout:
    fout.write(''.join(lines))
