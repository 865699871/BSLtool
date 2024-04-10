import os
import argparse

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-n", "--name_list_file")
    parser.add_argument("-r", "--genome_ref")
    parser.add_argument("-w", "--work_dir")
    parser.add_argument("-f", "--filter_identity")
    args = parser.parse_args()

    name_list_file = args.name_list_file
    genome_ref = args.genome_ref
    work_dir = args.work_dir
    filter_identity = float(args.filter_identity)

    script_path = os.path.split(os.path.realpath(__file__))[0]
    sd_path = script_path + '/stringdecomposer/bin/stringdecomposer'

    with open(name_list_file,'r') as nlf:
        while True:
            line = nlf.readline()[:-1]
            if not line:
                break
            items = line.split('_')
            chr = items[0]
            start = items[1]
            end = items[2]
            unit_len = int(items[3])
            # 先创建路径
            out_sd_dir = work_dir + '/' + line
            # samtools
            cmd = 'samtools faidx' + ' ' + genome_ref + ' ' + chr+':'+start+'-'+end + ' ' + '>' + out_sd_dir+'/'+line+'.fa'

            os.system(cmd)
            # sd
            overlap_threshold = 2 * unit_len
            if overlap_threshold > 500:
                # run SD
                cmd = 'python ' + sd_path + ' ' + out_sd_dir + '/' + line + '.fa' + ' ' + \
                      work_dir + '/' + line + '.fa' + \
                      ' -v ' + str(overlap_threshold)  + \
                      ' -o ' + out_sd_dir
            else:
                cmd = 'python ' + sd_path + ' ' + out_sd_dir + '/' + line + '.fa' + ' ' + \
                      work_dir + '/' + line + '.fa' + \
                      ' -o ' + out_sd_dir
            os.system(cmd)

            # get continue region
            sd_file = work_dir + '/' + line + '/final_decomposition.tsv'
            out_continue_region_file = work_dir + '/' + line + '/continue_region.txt'
            continuous_region = []
            region = []
            with open(sd_file, 'r') as sf:
                while True:
                    line = sf.readline()[:-1]
                    if not line:
                        break
                    items = line.split('\t')
                    start = int(items[2])
                    end = int(items[3])
                    identity = float(items[4])
                    if identity < filter_identity:
                        continue
                    region.append([start, end])

            if len(region) == 1:
                continuous_region.append(region[0])
            else:
                if len(region) != 0:
                    init_region = region[0]
                    for i in range(len(region) - 1):
                        if region[i + 1][0] - init_region[1] == 1:
                            init_region[1] = region[i + 1][1]
                        else:
                            continuous_region.append(init_region)
                            init_region = region[i + 1]
                    continuous_region.append(init_region)
            out_continue_region_file = open(out_continue_region_file, 'w')
            for i in continuous_region:
                out_continue_region_file.write(str(i[0]) + '\t' + str(i[1]) + '\n')
            out_continue_region_file.close()



if __name__ == '__main__':
    main()