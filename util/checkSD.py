import os
import argparse

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-n", "--name_list_file")
    parser.add_argument("-w", "--work_dir")
    args = parser.parse_args()

    name_list_file = args.name_list_file
    work_dir = args.work_dir

    with open(name_list_file,'r') as nlf:
        while True:
            line = nlf.readline()[:-1]
            if not line:
                break
            SD_dir = work_dir + '/' + line
            sd_file = SD_dir + '/' + 'final_decomposition.tsv'
            if not os.path.exists(sd_file):
                print(sd_file)




if __name__ == '__main__':
    main()