import argparse

from mysqlConnect.conn import conn
from greedyFinder.func import selecetdata
from greedyFinder.func import cover
from greedyFinder.func import hitcount
from greedyFinder.func import greedyfind
from resultView.func import pltview
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--subfamily", type=str, help="Target TE subfamily, default L1PA2",
                        default='L1PA2')
    parser.add_argument("-c", "--number_to_calculate", type=int, help="number of sgRNAs to calculate the combination, default 100",
                        default=100)
    parser.add_argument("-n", "--number", type=int, help="number of target sgRNAs, default 5", default=5)
    parser.add_argument("-t", "--threshold", type=float, help="threshold of coverage, default 0.8", default=0.8)
    parser.add_argument("-m", "--mode", type=float, help="greedy rank mode, 1: Strict mode (better for For evolutionary young TEs) 2: Lenient mode (better for ancient TEs)", default=1)
    args = parser.parse_args()
    n1 = args.number_to_calculate
    n2 = args.number
    threshold = args.threshold
    s = args.subfamily
    TE = '%' + s + ';%'
    mode = args.mode
    print("Loading TEs and potential sgRNAs...")
    df_TE, df_gRNA = selecetdata(TE, conn)
    # df_gRNA = pd.read_csv("df_gRNA.csv", sep=',', header=0, index_col=0, on_bad_lines='skip')
    # df_TE = pd.read_csv("df_TE.csv", sep=',', header=0, index_col=0, on_bad_lines='skip')
    print("Loading sgRNA hits...")
    df_sort_byCover, cover_dict = cover(df_TE, df_gRNA)
    df_sort_byCover = hitcount(df_sort_byCover, conn, n1)
    print("Searching for best sgRNAs combination...")
    coverage, set_gRNA, set_cover = greedyfind(df_sort_byCover, cover_dict, df_TE, n2, threshold, mode)
    print("Set of sgRNAs: ", set_gRNA)
    print("Coverage: ", coverage)
    pltview(coverage, set_gRNA, set_cover, df_sort_byCover)


