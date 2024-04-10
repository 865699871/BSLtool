data_dir = 'D:/Papaver_cen/buildsatellite_script/data/PB'

chr_number = 7

outfile = data_dir + '/' + 'TRF.bed'
outfile = open(outfile,'w')
for i in range(chr_number):
    chr = 'chr' + str(i + 1)
    data_file = data_dir + '/' + chr + '.fa.2.7.7.80.10.50.500.2.xls'
    with open(data_file,'r') as df:
        while True:
            line = df.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            start = items[0]
            end = items[1]
            RU = items[2]
            RN = items[3]
            unit = items[13]
            outfile.write(chr+'\t'+start+'\t'+end+'\t'+RU+'\t'+RN+'\t'+unit+'\n')
outfile.close()