import argparse

"""
Example: 
python3 merge_file.py -i LA_odd100_output_fst.txt -r LA_snplist.txt -q .1
"""

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--file', action='store', help='input file', required=True)
parser.add_argument('-r', '--ref', action='store', help='reference file', required=True)
parser.add_argument('-o', '--out', action='store', help='output file', default='output.txt')
parser.add_argument('-q', '--qval', action='store', type=float, default=.1, help='q value.')

args = parser.parse_args()
filename = args.file
reffile = args.ref
target = args.qval
outfile = args.out


with open(filename, 'r') as f, open(reffile, 'r') as reff, open(outfile, 'w') as fw:
	ls = f.readlines()
	rls = reff.readlines()
	headers = ls.pop(0).split()
	headers.insert(0, 'idx')

	fw.write('{0}\t{1}\t{2}\t{3}\n'.format('chr','pos', 'qval', 'fst'))
	
	try:
		assert len(ls) == len(rls)
	except:
		print('ERROR: file {0} ref file {1} length does not match!'.format(len(ls), len(rls)))

	cnt = 0
	for l, ref_l in zip(ls, rls):
		
		qval = float(l.split()[3])
		fst = float(l.split()[5])
		name = ref_l.strip('\n')

		if qval < target:
			cnt += 1
			fw.write('{0}\t{1}\t{2}\t{3}\n'.format(name.split('_')[0], name.split('_')[1], qval, fst))

	f.close()
	reff.close()

#print(cnt)
#print('Check \"{0}\" and \"{1}\".'.format(headers[3], headers[5]))
print('[{0}] records has q value < {1}.'.format(cnt, target))