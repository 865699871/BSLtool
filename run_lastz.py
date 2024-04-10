import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-n", "--name_list_file")
    parser.add_argument("-s", "--split_ref_dir")
    parser.add_argument("-w", "--work_dir")

    args = parser.parse_args()

    name_list_file = args.name_list_file
    split_ref_dir = args.split_ref_dir
    work_dir = args.work_dir

    filelist = os.listdir(split_ref_dir)

    with open(name_list_file,'r') as nlf:
        while True:
            line = nlf.readline()[:-1]
            if not line:
                break
            out_lastz_dir = work_dir + '/' + line
            if not os.path.exists(out_lastz_dir):
                os.mkdir(out_lastz_dir)
            target_fa = work_dir + '/' + line +'.fa'
            for i in filelist:
                ref = split_ref_dir + '/' + i
                cmd = 'lastz ' + ref + ' ' + target_fa + ' ' + \
                      '--format=general:score,name1,strand1,size1,start1,' \
                      'end1,name2,strand2,identity,' \
                      'length1,align1' + ' > ' + out_lastz_dir +'/' + i + '.xls'
                os.system(cmd)


if __name__ == '__main__':
    main()