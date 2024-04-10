import networkx as nx
import pandas as pd
import argparse
import os
def main():
    # 7.生成边文件
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-n", "--name_list_file")
    parser.add_argument("-w", "--workdir")
    parser.add_argument("-a", "--all_name_list_file")
    parser.add_argument("-s", "--split_ref_dir")

    parser.add_argument("-fi", "--filter_identity")
    parser.add_argument("-fa", "--filter_al")

    args = parser.parse_args()

    name_list_file = args.name_list_file
    workdir = args.workdir
    all_name_list_file = args.all_name_list_file
    split_ref_dir = args.split_ref_dir

    filter_identity = int(args.filter_identity)
    filter_al = int(args.filter_al)

    split_ref_list = os.listdir(split_ref_dir)

    print('start process')
    all_name_list = []
    with open(all_name_list_file, 'r') as nlf:
        while True:
            line = nlf.readline()[:-1]
            if not line:
                break
            all_name_list.append(line)

    name_info_table = {}
    count = 0
    for i in all_name_list:
        sd_continue_file = workdir + '/' + i + '/continue_region.txt'
        sd_continue = []
        with open(sd_continue_file,'r') as sf:
            while True:
                line = sf.readline()[:-1]
                if not line:
                    break
                items = line.split('\t')
                sd_continue.append([int(items[0]),int(items[1])])
        items = i.split('_')
        chr = items[0]
        start = int(items[1])
        end = int(items[2])
        unit_len = int(items[3])

        if chr not in name_info_table.keys():
            name_info_table[chr] = {}
            name_info_table[chr][i] = [chr, start, end, unit_len,sd_continue]
        else:
            name_info_table[chr][i] = [chr, start, end, unit_len,sd_continue]
        count += 1
    name_list = []
    with open(name_list_file, 'r') as nlf:
        while True:
            line = nlf.readline()[:-1]
            if not line:
                break
            name_list.append(line)

    edge_set = set()
    # 筛选条件，cov超过80%，且identity超过80
    count = 0
    print('start')
    for i in name_list:
        info = i.split('_')
        lastz_dir = workdir + '/' + i
        name_info = info
        name_lastz_set = set()

        for j in split_ref_list:
            lastz_chr_file = lastz_dir + '/' + j + '.xls'
            with open(lastz_chr_file,'r') as lcf:
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
                    if l_cov < int(name_info[-1]) * (filter_al/100):
                        continue
                    if l_identity < filter_identity:
                        continue
                    # 其他情况判断与node的关系
                    if l_chr not in name_info_table.keys():
                        continue

                    for k in name_info_table[l_chr].keys():
                        # 遍历区间看是否在其中
                        sd_continue = name_info_table[l_chr][k][4]
                        find = 0
                        name_info_k_start = int(k.split('_')[1])
                        for l in sd_continue:
                            n_start = name_info_k_start + int(l[0])
                            n_end = name_info_k_start+  int(l[1])
                            if n_start < l_start and l_end < n_end:
                                find = 1
                                name_lastz_set.add(k)
                                break
                        if find == 1:
                            break
        for j in name_lastz_set:
            if j == i:
                continue
            if (i+'$'+j in edge_set) or (j+'$'+i in edge_set):
                continue
            else:
                edge_set.add(i+'$'+j)
        count += 1
    index = name_list_file.split('/')[-1]
    edge_file = workdir+'/con_edge.'+index+'.txt'
    edge_file = open(edge_file,'w')
    for i in edge_set:
        edge_file.write(i+'\n')
    edge_file.close()




if __name__ == '__main__':
    main()