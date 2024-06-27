import matplotlib.pyplot as plt


def pltview(coverage, set_gRNA, set_cover, df_sort_byCover):
    # Plot components of all hits
    list_gRNA = list(set_gRNA)
    df_sel = df_sort_byCover[df_sort_byCover["grnaId"].isin(set_gRNA)]
    exonSum = df_sel["exonHit"].sum()
    ProSum = df_sel["promoterHit"].sum()
    total = len(set_cover) + exonSum + ProSum  # total hits
    final_hit = [len(set_cover), exonSum, ProSum]  # TEs Exon Promoter
    labels = ['TE subfamily', 'Exon', 'Promoter']
    colors = ['green', 'orange', 'yellow']
    sizes = [final_hit[0] / total * 100, final_hit[1] / total * 100, final_hit[2] / total * 100]
    expodes = (0.1, 0, 0)
    plt.pie(sizes, explode=expodes, labels=labels, shadow=True, colors=colors, autopct='%3.1f%%')
    plt.axis('equal')
    plt.savefig('Hits.png')
    plt.show()

    # Plot the total TE Coverage
    uncovered = 1 - coverage
    final_coverage = [coverage, uncovered]  # covered uncovered
    labels_1 = ['Covered TEs', 'Uncovered TEs']
    colors_1 = ['blue', 'grey']
    sizes_1 = [final_coverage[0], final_coverage[1]]
    expodes_1 = (0.1, 0)

    plt.pie(sizes_1, explode=expodes_1, labels=labels_1, shadow=True, colors=colors_1, autopct='%3.1f%%')
    plt.axis('equal')
    plt.savefig('Coverage.png')
    plt.show()
    return None