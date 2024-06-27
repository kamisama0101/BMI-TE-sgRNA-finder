import pandas as pd
from math import log, exp

def selecetdata(TE, conn):
    conn.ping(reconnect=True)
    # Get a cursor object that can execute the SQL command
    cursor = conn.cursor()  # By default, the result set returned after execution is displayed as a tuple
    # Get a cursor that can execute the SQL statement and return the result as a dictionary
    # cursor = conn.cursor(cursor=pymysql.cursors.DictCursor)

    # Extract TEs in the input subfamily from the database
    sql = "select chrom, chromStart, chromEnd, elementName, strand from dna_elements where elementName like '%s';" % TE
    # sql = "select * from dna_elements limit 5;"

    # The Command to execute an sql statement
    cursor.execute(sql)
    TE_locations = cursor.fetchall()
    # Initializes the list of gRNA locations
    gRNA_locations = []

    table = {'-': 'sgrna_negstrand', '+': 'sgrna_posstrand'}

    # Go through the sgRNA table in the Database to find all records of sgRNAs that hit the TE subfamilv
    for x in TE_locations:
        sql = "select chrom,chromStart,chromEnd,grnaId,strand from {table} force index(index_locations) where chrom='{chrom}' and chromStart>={start} and chromEnd<={end}" \
            .format(table=table.get(x[4]),
                    chrom=x[0],
                    start=x[1],
                    end=x[2])
        cursor.execute(sql)
        data = cursor.fetchall()
        gRNA_locations.extend(data)

    header_TE = ['chrom', 'chromStart', 'chromEnd', 'name', 'strand']
    header_gRNA = ['chrom', 'chromStart', 'chromEnd', 'grnaId', 'strand']
    df_TE = pd.DataFrame(list(TE_locations))
    df_gRNA = pd.DataFrame(list(gRNA_locations))
    df_TE.columns = header_TE[:len(df_TE.columns)]
    df_gRNA.columns = header_gRNA[:len(df_gRNA.columns)]

    # Write the information of TEs from the set subfamily and information of sgRNAs that target them
    # df_gRNA.to_csv('df_gRNA.csv')
    # df_TE.to_csv('df_TE.csv')

    # Unique the gRNAs according to their id
    list_gRNA = df_gRNA['grnaId'].unique()
    print("Number of possible sgRNAs: ", len(list_gRNA))
    print("Number of total hits on TEs: ", len(df_gRNA))

    # Close the cursor object
    cursor.close()

    # Closing the database connection
    conn.close()

    return df_TE, df_gRNA


def cover(df_TE, df_gRNA):
    cover_gRNA = df_gRNA.drop_duplicates(subset=['grnaId'], keep='first')
    cover_gRNA = cover_gRNA[["grnaId"]]
    cover_gRNA.index = range(len(cover_gRNA))
    cover_gRNA.insert(1, 'coverage', 0)
    cover_gRNA.insert(2, 'score', 0)
    cover_gRNA.insert(3, 'exonHit', 0)
    cover_gRNA.insert(4, 'promoterHit', 0)
    cover_dict = {}
    cover_dict[-1] = {-1}

    # Initialize the set of targets for each sgRNA
    for line in df_gRNA.itertuples(index=True):
        gid = line[4]
        cover_dict[gid] = set()

    # Use the soRNA No. as the index to get the set of targets for each sgRNA
    for line in df_TE.itertuples(index=True):
        tid = line[0]
        chrom = line[1]
        chromStart = line[2]
        chromEnd = line[3]
        df_t = df_gRNA[(df_gRNA['chrom'] == chrom) & (df_gRNA['chromStart'] >= chromStart) & (df_gRNA['chromEnd'] <= chromEnd)]
        for x in zip(df_t['grnaId']):
            gid = x[0]
            cover_dict[gid].add(tid)
    # Get the rank of each sgRNA according to the TE Coverage from most to least
    for line in df_gRNA.itertuples(index=True):
        gid = line[4]
        cover_dict[gid].discard(-1)  # Deletes the initialized value
        cover_gRNA.loc[cover_gRNA['grnaId'] == gid, 'coverage'] = len(cover_dict[gid])
    df_sort_byCover = cover_gRNA.sort_values(by="coverage", ascending=False)

    return df_sort_byCover, cover_dict


def hitcount(df_sort_byCover, conn, n1):

    # Connect To Database
    conn.ping(reconnect=True)
    cursor = conn.cursor()

    # Traverse the genome annotation file to get the number of exons or enhancers targeted by each sgRNA
    sql1 = "select * from {table} force index(index_sgrnaId) where grnaId={gid}"
    sql2 = "select elementName from dna_elements force index(index_locations) where chrom='{chrom}' and chromStart<={left} and chromEnd>={right} and strand='{strand}'"
    t = 0

    for line in df_sort_byCover.itertuples(index=True):
        result_set = set()
        gid = line[1]
        sql = sql1.format(
            table='sgrna_negstrand',
            gid=gid
        )
        cursor.execute(sql)
        dataneg = cursor.fetchall()
        sql = sql1.format(
            table='sgrna_posstrand',
            gid=gid
        )
        cursor.execute(sql)
        datapos = cursor.fetchall()
        data = datapos + dataneg

        for x in data:
            chrom = x[0]
            left = x[1]
            right = x[2]
            strand = x[5]
            sql = sql2.format(
                chrom=chrom,
                left=left,
                right=right,
                strand=strand
            )
            cursor.execute(sql)
            hit = cursor.fetchall()
            for p in hit:
                result_set.add(p[0])

        for x in data:
            chrom = x[0]
            left = x[1]
            right = x[2]
            strand = x[5]
            sql = sql2.format(
                chrom=chrom,
                left=left,
                right=right,
                strand=strand
            )
            cursor.execute(sql)
            hit = cursor.fetchall()
            for p in hit:
                result_set.add(p[0])

        for i in result_set:
            if 'promoter-TSS' in i:
                df_sort_byCover.loc[df_sort_byCover['grnaId'] == gid, 'promoterHit'] += 1
            if 'exon' in i:
                df_sort_byCover.loc[df_sort_byCover['grnaId'] == gid, 'exonHit'] += 1
        t += 1
        if t == n1:
            break
    cursor.close()
    conn.close()

    # Get 100 sgRNAs with top coverages of TEs
    df_sort_byCover = df_sort_byCover.iloc[0:n1, ]
    # df_sort_byCover.to_csv('df_sort_byCover.csv')

    return df_sort_byCover


def greedyfind(df_sort_byCover, cover_dict, df_TE, n2, threshold, mode):
    # The scoring mechanism
    if mode == 1:
        df_sort_byCover["score"] = df_sort_byCover["coverage"] / (
                2 * df_sort_byCover["exonHit"] + df_sort_byCover["promoterHit"] + 1)
    elif mode == 2:
        df_sort_byCover["score"] = df_sort_byCover["coverage"] / (
                log(df_sort_byCover["exonHit"] + 10, 10) + log(df_sort_byCover["promoterHit"] + 10, 10))
    # Get the rank of each sgRNA by considering exons and TSS-promoters
    df_sort_byCover = df_sort_byCover.sort_values(by="score", ascending=False)
    set_cover = {-1}
    set_gRNA = {-1}
    coverage = 0

    # Loop until there are 10 sgRNAs in the result set or coverage greater than 0.9
    while len(set_gRNA) <= n2 and coverage <= threshold:
        nextId = -1
        for line in df_sort_byCover.itertuples(index=True):
            gid = line[1]
            if gid not in set_gRNA:
                nextId = gid
                break
        set_gRNA.add(nextId)
        set_cover = set_cover.union(cover_dict[nextId])
        for line in df_sort_byCover.itertuples(index=True):
            gid = line[1]
            x = len(set_cover.union(cover_dict[gid])) - len(set_cover)
            y = 1
            if mode == 1:
                y = 2 * df_sort_byCover.loc[df_sort_byCover["grnaId"] == gid, "exonHit"] + df_sort_byCover.loc[
                    df_sort_byCover["grnaId"] == gid, "promoterHit"]
            elif mode == 2:
                y = log(df_sort_byCover.loc[df_sort_byCover["grnaId"] == gid, "exonHit"] + 10, 10) + log(df_sort_byCover.loc[
                    df_sort_byCover["grnaId"] == gid, "promoterHit"] + 10, 10)
            df_sort_byCover.loc[df_sort_byCover["grnaId"] == gid, "score"] = x / y
        df_sort_byCover = df_sort_byCover.sort_values(by="score", ascending=False)
        coverage = len(set_cover) / len(df_TE)

    set_gRNA.discard(-1)
    set_cover.discard(-1)
    return coverage, set_gRNA, set_cover