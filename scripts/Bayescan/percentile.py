import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--file', action='store', help='input file', required=True)
parser.add_argument('-o', '--output', action='store', help='output file', default='output.txt')
parser.add_argument('-p', '--percentile', action='store', type=float, default=99, help='indicated percentile. default 99.')

dat = []
args = parser.parse_args()
filename = args.file
print('Open file ' + filename + '...')
percentile = args.percentile
outfile = args.output

f = open(filename, 'r')
Header = True
header = ''
cnt = 0
for l in f:
    if Header:
        header = l
        Header = False
        continue
    _d = l.split()[3]
    dat.append(float(_d))
    cnt += 1
f.close()

print('Read ' + str(cnt) + ' records.')
target = np.percentile(dat, percentile)

print('Percentile target for percentile {0} is {1}'.format(percentile, target))
print('Max value is: ' + str(max(dat)))
print('Now extract samples > target...')

f = open(filename, 'r')

with open(outfile, 'w') as fw:
    fw.write(header)
    tcnt = 0
    Header = True
    for l in f:
        if Header:
            Header = False
            continue
        _d = float(l.split()[3])
        if _d > target:
            tcnt += 1
            fw.write(l)

f.close()

print('Total found {0} records > {1}.'.format(tcnt, target))
