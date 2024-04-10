import argparse
import os
import networkx as nx
from joblib import Parallel, delayed
import networkx.algorithms.community as nx_comm

def processTRFresult(trf_result_file,
                     filter_repeat_number,
                     filter_repeat_unit_length,
                     overlap_ratio):
    init_trf = {}
    with open(trf_result_file, 'r') as trf_f:
        while True:
            line = trf_f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            chr = items[0]
            start = int(items[1])
            end = int(items[2])
            repeat_unit_length = int(items[3])
            repeat_number = float(items[4])
            repeat_unit = items[5]
            # filter repeat number < 100 and repeat unit length < 100
            if (repeat_unit_length > filter_repeat_unit_length) & (repeat_number > filter_repeat_number):
                if chr not in init_trf.keys():
                    init_trf[chr] = [[start, end, repeat_unit_length, repeat_unit]]
                else:
                    init_trf[chr].append([start, end, repeat_unit_length, repeat_unit])

    # filter overlap > overlap ratio
    filter_trf = {}
    for i in init_trf.keys():
        chr = i
        filter_trf[chr] = []
        sorted_region = sorted(init_trf[i], key=lambda x: x[0])
        start_one = sorted_region[0]
        # if next region overlap with current region larger than overlap ratio,
        # save small unit length one
        for j in range(len(sorted_region) - 1):
            index = j + 1
            if start_one[1] >= sorted_region[index][0]:
                overlap = start_one[1] - sorted_region[index][0]

                if overlap / (sorted_region[index][1] - sorted_region[index][0]) > overlap_ratio:
                    if start_one[2] > sorted_region[index][2]:
                        start_one = sorted_region[index]
                    else:
                        start_one = start_one
                else:
                    filter_trf[chr].append(start_one)
                    start_one = sorted_region[index]
            else:
                filter_trf[chr].append(start_one)
                start_one = sorted_region[index]
        filter_trf[chr].append(start_one)

    TRF_result = []
    for i in filter_trf.keys():
        for j in filter_trf[i]:
            TRF_result.append([i, j[0], j[1], j[2], j[3]])

    return TRF_result


def splitTRF(TRF_result, threads, outdir):
    one_file_region_number = int(len(TRF_result) / threads) + 1
    count = 0
    file_index = 1
    outfile = outdir + '/split.' + str(file_index) + '.txt'
    file_list = [outfile]
    outfile = open(outfile, 'w')
    all_regions = []
    for i in TRF_result:
        name = i[0] + '_' + str(i[1]) + '_' + str(i[2]) + '_' + str(i[3])
        all_regions.append(name)
        outfafile = outdir + '/' + name + '.fa'
        outfafile = open(outfafile, 'w')
        outfafile.write('>' + name + '\n')
        outfafile.write(i[4] + '\n')
        outfafile.close()
        if not os.path.exists(outdir + '/' + name):
            os.mkdir(outdir + '/' + name)

        if count == one_file_region_number:
            outfile.close()
            file_index += 1
            outfile = outdir + '/split.' + str(file_index) + '.txt'
            file_list.append(outfile)
            outfile = open(outfile, 'w')
            count = 0
            outfile.write(name + '\n')
        else:
            outfile.write(name + '\n')
            count += 1
    outfile.close()
    return file_list, all_regions


def run_SD(region_file, SD_script, genome_ref, filter_identity, out_region_dir):
    cmd = 'python ' + \
          SD_script + ' ' + \
          '-n' + ' ' + region_file + ' ' + \
          '-r' + ' ' + genome_ref + ' ' + \
          '-w' + ' ' + out_region_dir + ' ' + '-f' + ' ' + str(filter_identity)
    os.system(cmd)
    print(cmd)


def run_lastz(region_file, lastz_script, split_ref_dir, out_region_dir):
    cmd = 'python ' + \
          lastz_script + ' ' + \
          '-n' + ' ' + region_file + ' ' + \
          '-s' + ' ' + split_ref_dir + ' ' + \
          '-w' + ' ' + out_region_dir
    os.system(cmd)
    print(cmd)


def run_buildEdges(region_file, buildEdges_script, out_region_dir, all_regions_filename, split_ref_dir, filter_identity,
                   filter_al):
    cmd = 'python ' + \
          buildEdges_script + ' ' + \
          '-n' + ' ' + region_file + ' ' + \
          '-w' + ' ' + out_region_dir + ' ' + \
          '-a' + ' ' + all_regions_filename + ' ' + \
          '-s' + ' ' + split_ref_dir + ' ' + \
          '-fi' + ' ' + str(filter_identity) + ' ' + \
          '-fa' + ' ' + str(filter_al)
    os.system(cmd)
    print(cmd)

def satellite_lastz(i,out_satellites_dir,out_split_region_dir,split_ref_dir):
    satellite_dir = out_satellites_dir + '/' + i
    if not os.path.exists(satellite_dir):
        os.mkdir(satellite_dir)

    templet_fa = out_split_region_dir + '/' + i + '.fa'
    # run lastz
    split_ref_filelist = os.listdir(split_ref_dir)
    for j in split_ref_filelist:
        ref = split_ref_dir + '/' + j
        cmd = 'lastz ' + ref + ' ' + templet_fa + ' ' + \
              '--format=general:score,name1,strand1,size1,start1,' \
              'end1,name2,strand2,identity,' \
              'length1,align1' + ' > ' + satellite_dir + '/' + j + '.xls'
        os.system(cmd)


def main():
    parser = argparse.ArgumentParser(description="Building satellite library")
    parser.add_argument("-i", "--trf_result_file", help="TRF result bed, required", required=True)
    parser.add_argument("-g", "--genome_ref", help="genome, required", required=True)
    parser.add_argument("-o", "--outdir", help="HiCAT reads output path default is ./BSL_out",
                        default='./BSL_out',
                        required=False)

    parser.add_argument("-th", "--threads", help="The number of threads, default is 1", type=int, default=1,
                        required=False)

    parser.add_argument("-rn", "--filter_repeat_number", help="filter repeat number, default is 100", type=int,
                        default=100,
                        required=False)
    parser.add_argument("-ru", "--filter_repeat_unit_length", help="filter repeat unit length, default is 100",
                        type=int,
                        default=100,
                        required=False)
    parser.add_argument("-f", "--filter_identity", help="filter identity, default is 80%%",
                        type=int,
                        default=80,
                        required=False)
    parser.add_argument("-fa", "--filter_al", help="filter alignment length in lastz, default is 80%%",
                        type=int,
                        default=80,
                        required=False)
    parser.add_argument("-m", "--merge_interval", help="merge interval for satellite regions, default is 500K",
                        type=int,
                        default=500000,
                        required=False)
    parser.add_argument("-or", "--overlap_ratio", help="filter overlap ratio in TRF, default is 0.6",
                        type=float,
                        default=0.6,
                        required=False)

    args = parser.parse_args()

    trf_result_file = args.trf_result_file
    genome_ref = args.genome_ref
    outdir = args.outdir

    threads = args.threads

    filter_repeat_number = args.filter_repeat_number
    filter_repeat_unit_length = args.filter_repeat_unit_length
    filter_identity = args.filter_identity
    filter_al = args.filter_al
    merge_interval = args.merge_interval
    overlap_ratio = args.overlap_ratio

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    script_path = os.path.split(os.path.realpath(__file__))[0]
    SD_script = script_path + '/run_SD.py'
    lastz_script = script_path + '/run_lastz.py'
    buildEdges_script = script_path + '/buildEdge.py'
    TRF_result = processTRFresult(trf_result_file,
                                  filter_repeat_number,
                                  filter_repeat_unit_length,
                                  overlap_ratio)

    if len(TRF_result) != 0:
        out_split_region_dir = outdir + '/split_region'
        if not os.path.exists(out_split_region_dir):
            os.mkdir(out_split_region_dir)
        file_list, all_regions = splitTRF(TRF_result, threads, out_split_region_dir)
        # string decomposer
        Parallel(n_jobs=threads)(delayed(run_SD)(i,
                                                 SD_script,
                                                 genome_ref,
                                                 filter_identity,
                                                 out_split_region_dir) for i in file_list)
        # split ref and lastz
        chr_list = set()
        for i in TRF_result:
            chr_list.add(i[0])
        split_ref_dir = outdir + '/split_ref'
        if not os.path.exists(split_ref_dir):
            os.mkdir(split_ref_dir)
        for i in chr_list:
            cmd = 'samtools faidx' + ' ' + genome_ref + ' ' + i + ' ' + \
                  '>' + split_ref_dir + '/' + i + '.fa'
            os.system(cmd)
        Parallel(n_jobs=threads)(
            delayed(run_lastz)(i,
                               lastz_script,
                               split_ref_dir,
                               out_split_region_dir) for i in file_list)

        # get edges
        all_regions_filename = out_split_region_dir + '/all_regions.txt'
        all_regions_file = open(all_regions_filename, 'w')
        for i in all_regions:
            all_regions_file.write(i + '\n')
        all_regions_file.close()

        Parallel(n_jobs=threads)(delayed(run_buildEdges)(i,
                                                         buildEdges_script,
                                                         out_split_region_dir,
                                                         all_regions_filename,
                                                         split_ref_dir,
                                                         filter_identity,
                                                         filter_al) for i in file_list)

        # merge edges build network
        edge_files = []
        all_files = os.listdir(out_split_region_dir)
        for i in all_files:
            if i.startswith('con_edge'):
                edge_files.append(i)
        edge_set = set()
        for i in edge_files:
            edge_file = out_split_region_dir + '/' + i
            with open(edge_file, 'r') as ef:
                while True:
                    line = ef.readline()[:-1]
                    if not line:
                        break
                    items = line.split('$')
                    if line in edge_set:
                        continue
                    if items[1] + '$' + items[0] in edge_set:
                        continue
                    edge_set.add(line)
        outedgefile = out_split_region_dir + '/edges.xls'
        outedgefile = open(outedgefile, 'w')
        G = nx.Graph()
        for i in edge_set:
            items = i.split('$')
            G.add_edge(items[0], items[1])
            outedgefile.write(i + '\n')
        outedgefile.close()
        community = nx_comm.louvain_communities(G, seed=1)

        community_set = set()
        for i in community:
            for j in i:
                community_set.add(j)

        out_community_file = out_split_region_dir + '/community.xls'
        out_community_file = open(out_community_file, 'w')
        for i in list(community):
            # 找unit_len最多的中度大的节点作为community name
            region_number = len(i)
            unit_table = {}
            for j in i:
                items = j.split('_')
                unit_len = int(items[-1])
                if unit_len not in unit_table.keys():
                    unit_table[unit_len] = 1
                else:
                    unit_table[unit_len] += 1
            max_degree = -1
            max_degree_name = ''
            max_unit_len = -1
            max_unit = -1
            for j in unit_table.keys():
                if max_unit_len < unit_table[j]:
                    max_unit_len = unit_table[j]
                    max_unit = j
            for j in i:
                items = j.split('_')
                unit_len = int(items[-1])

                if unit_len == max_unit:
                    if G.degree(j) > max_degree:
                        max_degree = G.degree(j)
                        max_degree_name = j
            out_community_file.write(max_degree_name)
            out_community_file.write('\t' + str(max_unit))
            out_community_file.write('\t' + str(region_number))
            for j in i:
                out_community_file.write('\t' + j)
            out_community_file.write('\n')

        for i in all_regions:
            if i in community_set:
                continue
            out_community_file.write(i)
            info = i.split('_')
            unit_len = int(info[-1])
            region_number = 1
            out_community_file.write('\t' + str(unit_len))
            out_community_file.write('\t' + str(region_number))
            out_community_file.write('\t' + i)
            out_community_file.write('\n')

        out_community_file.close()

        # build satellite library rank based on AT and lastz to get region
        satelite_communities = []
        with open(out_split_region_dir + '/community.xls','r') as cf:
            while True:
                line = cf.readline()[:-1]
                if not line:
                    break
                items = line.split('\t')[0]
                satelite_communities.append(items)
        out_satellites_dir = outdir + '/satellites'
        if not os.path.exists(out_satellites_dir):
            os.mkdir(out_satellites_dir)
        print('satellite lastz')
        Parallel(n_jobs=threads)(delayed(satellite_lastz)(i,
                                                          out_satellites_dir,
                                                          out_split_region_dir,
                                                          split_ref_dir) for i in satelite_communities)
        split_ref_filelist = os.listdir(split_ref_dir)
        satellite_content = []
        for i in satelite_communities:
            templet_fa = out_split_region_dir + '/' + i + '.fa'
            satellite_dir = out_satellites_dir + '/' + i
            templet_seq = ''
            with open(templet_fa,'r') as tf:
                tf.readline()
                templet_seq = tf.readline()[:-1]
            AT = 0
            for j in templet_seq:
                if j == 'A' or j == 'T':
                    AT += 1
            AT_content = AT / len(templet_seq)
            # genome size
            satellite_regions = []
            satellite_regions_merge = []
            genome_size = 0
            satellite_chr = set()
            for j in split_ref_filelist:
                satellite_lastz_file = satellite_dir + '/' + j + '.xls'
                chr_regions = []
                with open(satellite_lastz_file, 'r') as lcf:
                    while True:
                        line = lcf.readline()[:-1]
                        if not line:
                            break
                        if line.startswith('#'):
                            continue
                        items = line.split('\t')
                        l_chr = items[1]
                        l_start = int(items[4])
                        l_end = int(items[5])
                        l_cov = int(items[8].split('/')[-1])
                        l_identity = float(items[9][:-1])
                        if l_cov < len(templet_seq) * (filter_al / 100):
                            continue
                        if l_identity < filter_identity:
                            continue
                        genome_size += (l_end - l_start)
                        satellite_regions.append([l_chr,l_start,l_end])
                        chr_regions.append([l_chr,l_start,l_end])
                        if l_chr not in satellite_chr:
                            satellite_chr.add(l_chr)
                # merge region
                if len(chr_regions) != 0:
                    sorted_chr_regions = sorted(chr_regions,key=lambda x:x[1])

                    init_start = sorted_chr_regions[0]
                    for j in range(len(sorted_chr_regions) - 1):
                        if sorted_chr_regions[j + 1][1] - init_start[2] < merge_interval:
                            init_start[2] = sorted_chr_regions[j + 1][2]
                        else:
                            satellite_regions_merge.append([init_start[0], init_start[1], init_start[2]])
                            init_start = sorted_chr_regions[j + 1]
                    satellite_regions_merge.append([init_start[0], init_start[1], init_start[2]])

            # index (rank by AT), represented region,
            # chromosome, repeat number, AT content, genome size, templet seq
            satellite_content.append([i,
                                      satellite_chr,
                                      len(satellite_regions),
                                      AT_content,
                                      genome_size,
                                      AT_content * genome_size,
                                      templet_seq,
                                      satellite_regions,
                                      satellite_regions_merge])

        sorted_satellites = sorted(satellite_content,key=lambda x:x[5],reverse=True)
        outfile = out_satellites_dir + '/summary_satellite.xls'
        outfile = open(outfile,'w')
        count = 1
        for i in sorted_satellites:
            represented_region = i[0]
            satellite_chr = i[1]
            repeat_number = i[2]
            AT_content = i[3]
            genome_size = i[4]
            templet_seq = i[6]

            outfile.write(str(count) +'\t'+represented_region+'\t')
            chr_str = ''
            for j in satellite_chr:
                chr_str += j+','
            chr_str = chr_str[:-1]
            outfile.write(chr_str+'\t'+str(repeat_number)+'\t'+str(AT_content)+'\t'+str(genome_size)+'\t'+templet_seq+'\n')

            satellite_regions = i[7]
            out_satellite_region_file = out_satellites_dir + '/' + str(count) + '_' + represented_region+'.region.bed'
            out_satellite_region_file = open(out_satellite_region_file,'w')
            for j in satellite_regions:
                out_satellite_region_file.write(j[0]+'\t'+str(j[1])+'\t'+str(j[2])+'\n')
            out_satellite_region_file.close()
            satellite_regions_merge = i[8]
            out_satellite_region_merge_file = out_satellites_dir + '/' + str(count) + '_' + represented_region+'.mergeregion'+str(merge_interval)+'.bed'
            out_satellite_region_merge_file = open(out_satellite_region_merge_file, 'w')
            for j in satellite_regions_merge:
                out_satellite_region_merge_file.write(j[0]+'\t'+str(j[1])+'\t'+str(j[2]) +'\t' +str(j[2] - j[1])+'\n')
            out_satellite_region_merge_file.close()
            count += 1
        outfile.close()

    else:
        print('no TRF result')


if __name__ == '__main__':
    main()
