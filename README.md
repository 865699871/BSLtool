# BSLtool: Build satellite library
This pipline is designed for building satellite library based on TRF method. It includes three steps (TRF results are input). First, BSLtool performed StringDecomposer based on the consensus repeat sequence to remove the noise in each TRF region. Then, using Lastz to build connection between each repeat region. Finally, BSLtool performed community detection method for TR network to build satellite library.
![img.png](img.png)
## Dependencies
Python 3.9.16

Packages  | Version used in Research|
--------- | --------|
joblib  | 1.1.0 |
lastz  | 1.04.22 |
samtools  | 1.6 |
networkx  | 2.7.1 |

StringDecomposer (https://github.com/ablab/stringdecomposer) with version 1.1.2.   
Development environment: Linux  
Development tool: Pycharm
## Installation

#### Source code (g++ version 5.3.1 or higher for stringdecomposer)

```
#install
conda install -y --file requirements.txt
cd ./stringdecomposer && make
```

## Usage
```Bash
usage: buildSatelliteLibrary.py [-h] -i TRF_RESULT_FILE -g GENOME_REF [-o OUTDIR] [-th THREADS] [-rn FILTER_REPEAT_NUMBER] [-ru FILTER_REPEAT_UNIT_LENGTH] [-f FILTER_IDENTITY]
                                [-fa FILTER_AL] [-m MERGE_INTERVAL] [-or OVERLAP_RATIO]

Building satellite library

optional arguments:
  -h, --help            show this help message and exit
  -i TRF_RESULT_FILE, --trf_result_file TRF_RESULT_FILE
                        TRF result bed, required
  -g GENOME_REF, --genome_ref GENOME_REF
                        genome, required
  -o OUTDIR, --outdir OUTDIR
                        HiCAT reads output path default is ./BSL_out
  -th THREADS, --threads THREADS
                        The number of threads, default is 1
  -rn FILTER_REPEAT_NUMBER, --filter_repeat_number FILTER_REPEAT_NUMBER
                        filter repeat number, default is 100
  -ru FILTER_REPEAT_UNIT_LENGTH, --filter_repeat_unit_length FILTER_REPEAT_UNIT_LENGTH
                        filter repeat unit length, default is 100
  -f FILTER_IDENTITY, --filter_identity FILTER_IDENTITY
                        filter identity, default is 80%
  -fa FILTER_AL, --filter_al FILTER_AL
                        filter alignment length in lastz, default is 80%
  -m MERGE_INTERVAL, --merge_interval MERGE_INTERVAL
                        merge interval for satellite regions, default is 500K
  -or OVERLAP_RATIO, --overlap_ratio OVERLAP_RATIO
                        filter overlap ratio in TRF, default is 0.6

```
Inputs

TRF_RESULT_FILE e.g. TRF.bed (chr start end unit_length consensus with \t separator) 
```Bash
chr1    1       13378   7       1911.1  CCCTGAA
chr1    24573   24924   7       50.4    TGAACCC
chr1    95610   95690   1       81.0    A
chr1    95610   95690   31      2.5     AAAAAAAGAAAAAAAAATAAAAAAAAAAAAA
chr1    95610   95690   34      2.4     AAAAAAAGAAAAAAAAATCAAAAAAAAAAAAAGA
chr1    153388  153426  13      3.0     AGTTCGCGGACTG
chr1    153485  153523  13      3.0     GTTCGCGGACTGG
chr1    159345  159383  13      3.0     GTTCGCGGACTGG
chr1    182684  182709  13      2.0     CCAGTCCGCGAAC
chr1    182769  182807  13      3.0     CAGTCCGCGAACT
chr1    186645  194018  115     64.4    CATCACCCTACGCCTCTGCCAAGCAATAGAGCAAACTGGGAATCACCCCACTTGTCATATTGCTAAGATTTGGATTACATGTTCACGGGTTGGATTGACGGCAGATGAGACTAGG
chr1    199043  199143  33      3.1     AACAAAAGACGGATTCGGCGATCAAGGCTTTGG
chr1    205863  205913  26      2.0     CACCCCTCACAACGATCGGTTCCCTA
chr1    205972  206073  27      3.8     CCCCATACAAGGATCGGTTCCCTTTCA
chr1    205972  206075  53      2.0     CCCCACACAAGGATCGATTCCCTTTCACCCCATACAAGGATCGGTTCCCATCA
...

```
GENOME_REF
```Bash
>chr1
CCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAAC
CCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACC
CAGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCC
TGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCT
GAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTG
AACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGA
ACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAA
......
>chr2
CTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCC
TGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCT
GAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTG
AACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGA
ACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAA
CCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAACCCTGAAC
......
```
Other Parameters(default)
```Bash
Regions with unit length lower than FILTER_REPEAT_UNIT_LENGTH(100) and repeat number lower than FILTER_REPEAT_NUMBER(100) will be removed.
Regions overlap with others larger than OVERLAP_RATIO(0.6) will be removed and saved the region with small unit.
Each region will be denoised by StringDecomposer. The filter identity is FILTER_IDENTITY(80%).
Regions with identity large than FILTER_IDENTITY(80%) and alignment length larger than FILTER_AL(80%) will be considered as same satellite.
BSLtool merged the closed regions with distance smaller than MERGE_INTERVAL(500000)
```

### Output
```Bash
The output is ./BSL_out/satellite
summary_satellite.xls (index default_name chrs repeat_number AT_percentage genomic_size unit_sequence)
1       chr4_5124774_5258915_196        chr2,chr3,chr6,chr1,chr4,chr7,chr5      406849  0.6020408163265306      79482756        CAAATAGAGACCAAAGGAGTTGCAAAAGGATATTGAGAGGTATAAGGTATCGAAGAGC
2       chr3_341478431_341708947_165    chr2,chr3,chr6,chr1,chr4,chr7,chr5      255554  0.6121212121212121      41974357        TATCTGGGTGGAAACATCTAACTAATAAAAGAGACCGCGGGCACCAATAGAAATTATT
3       chr3_353718464_353764854_371    chr2,chr3,chr6,chr1,chr4,chr7,chr5      83809   0.6091644204851752      31640855        ACATGTTGATCTGTAGATCCACATGGTCACCAGATTTTAGCCAAATTCCATCGTTTAG
4       chr1_51500271_51547455_162      chr2,chr3,chr6,chr1,chr4,chr7,chr5      93698   0.61875 14691166        AAATATATCAGAATTGAACATCACAGGTTCAAAGATTACTAGAGTTAGACCTGAGTTTCTCGTTCGGATGAACT
5       chr4_34804587_34821359_112      chr2,chr3,chr6,chr1,chr4,chr7,chr5      80052   0.6875  8611586 TTACTTTACTTATTTCCGATAAGGGTAGTATGGTCATTTTCACAGTCTCCATAATATTTATTGCACTTATTTCCAATAAGGG
6       chr5_240409096_240502546_194    chr2,chr3,chr6,chr1,chr4,chr7,chr5      51128   0.5958549222797928      9866902 GATCAATTTTGAGATTCCCAAGACCGTAGGCTTGTGCTTTCATGCCCATGCCCAAGATTCAATTTG
7       chr6_82151054_82215272_238      chr2,chr3,chr6,chr1,chr4,chr7,chr5      23525   0.7276595744680852      5409299 CATAGCATTGGTATCTACCAGCTAGAAGATTGTTTAAACCAACTGAAAGTAACTTAATTCCTTATA
8       chr1_173987628_174011000_169    chr1,chr5       11354   0.6863905325443787      1901624 AGTATAGTCAATTGCTTTGTTTTTATAGAAATTTTATTATATAACTGTGGAAGAGTCAAAAAATAAATATTCTTACAGAGGCCACCACGA
9       chr3_174320291_174340918_151    chr2,chr3,chr6,chr1,chr4,chr7,chr5      11514   0.5592105263157895      1717470 TTCCTCAACACATCCTTCCTCTTAGTGTTGAGACGTGTAAATTCCTCAACACTCACTACTGTGTTG
10      chr5_261623666_261637016_122    chr2,chr3,chr6,chr1,chr4,chr7,chr5      13982   0.4672131147540984      1672920 CGCTATTTTTGGCATGTGCCAAACGGTATAGCGGCTTAGAAATTGGCTAGTGGTGCGTGTTGGCAT
11      chr1_16467589_16597481_152      chr2,chr3,chr6,chr1,chr4,chr7,chr5      9098    0.5789473684210527      1283863 AGAACGATGAATCTATCCTTTTACAAGTAAGATTCCATTGGAATATTCTCCGAAACCTTAAACACA
12      chr7_173082461_173125513_418    chr2,chr3,chr1,chr4,chr7,chr5   3144    0.5047846889952153      1311172 ACTTGTAGAAATGCGATTAATACATATATAGGGAAGGGGTGTTGGCTTGAGCTTGAAAAATAAAGTTGGTCCC
......
*.region.bed:
Final regions for each satellite.
*.mergeregionMERGE_INTERVAL.bed: 
Final regions merged with distance smaller than MERGE_INTERVAL(500000). For visualization or HOR analysis(e.g. HiCAT)
```


## Contact
If you have any questions, please feel free to contact: gaoxian15002970749@163.com, xfyang@xjtu.edu.cn, kaiye@xjtu.edu.cn





