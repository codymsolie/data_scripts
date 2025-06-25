# Data_Generator.sage
#
# drives the logic of the program. 
# calls upon Common, EKR, Graph_Properties,
# and Intersection_Density classes to 
# collect relevant data, and saves it 
# to Neon. 

import os
from psycopg2 import pool
from dotenv import load_dotenv

def db_save(data, eigenvalues):
    load_dotenv()
    connection_string = os.getenv('DATABASE_URL')
    connection_pool = pool.SimpleConnectionPool(
        1,  # Minimum number of connections in the pool
        10,  # Maximum number of connections in the pool
        connection_string
    )
    if connection_pool:
        print("neondb connection successful.")
    conn = connection_pool.getconn()
    cur = conn.cursor()

    cur.execute("SELECT group_id FROM ekr_data WHERE \
                degree=%s AND gap_id=%s;",
                (data[1], data[2]))
    result = cur.fetchone()

    if (result):
        # update existing group
        # i dont feel the need to update eigenvalues. the computation will never change.
        print("updating data")
        cur.execute("UPDATE ekr_data SET size=%s, struc_desc=%s, int_dens_hi=%s, \
                     int_dens_lo=%s,  int_dens=%s, transitivity=%s, min_trans=%s,\
                     is_union=%s, is_join=%s, is_cmp=%s, is_pm_join=%s, is_cograph=%s, \
                     ekr=%s, is_abelian=%s, is_nilpotent=%s, is_primitive=%s\
                     WHERE degree=%s AND gap_id=%s;",
                    (data[3], data[4], data[5], data[6], data[7], data[8], data[9], data[10],
                     data[11], data[12], data[13], data[14], data[15], data[16], data[17], 
                     data[18], data[1], data[2]))
        conn.commit()

    else:
        # new group
        print("new data")
        cur.execute("INSERT INTO ekr_data \
                    (name, degree, gap_id, size, struc_desc, int_dens_hi, int_dens_lo, int_dens, \
                    transitivity, min_trans, is_union, is_join, is_cmp, is_pm_join, is_cograph,  \
                    ekr, is_abelian, is_nilpotent, is_primitive) \
                    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);",
                    data)
        conn.commit()

        # get previous group id to link eigenvalues
        cur.execute("SELECT group_id FROM ekr_data WHERE degree=%s AND gap_id=%s;",
                    (data[1], data[2]))
        gid=cur.fetchone()

        for eigenvalue in eigenvalues:
            cur.execute("INSERT INTO eigenvalues (group_id, eigenvalue, multiplicity) \
                        VALUES (%s, %s, %s)",
                        (gid, eigenvalue[0], eigenvalue[1]))
            conn.commit()

    # cleanup cleanup everybody do their share
    cur.close()
    connection_pool.putconn(conn)
    connection_pool.closeall()
    return

class Data_Generator:
    def __init__(self, groups):
        groups_left = len(groups)
        for G in groups:
            print("groups left: ", groups_left)
            print(G)
            common = Common(G)
            graph = Graph_Properties(common)
            ekr = EKR(common)
            int_dens = Intersection_Density(common, graph, ekr)
            groups_left -= 1

            data = (
                str(G),
                int(common.degree),
                int(gap.TransitiveIdentification(G)),
                int(common.order),
                str(G.structure_description()),
                float(int_dens.upper_bound),
                float(int_dens.lower_bound),
                float(int_dens.exact_value),
                int(common.transitivity),
                bool(common.minimally_transitive),
                bool(graph.is_union),                       #
                bool(graph.is_join),
                bool(graph.is_cmp),
                bool(graph.is_pm_join),                     #
                bool(graph.is_cograph),                     #
                int(ekr.has_ekr),
                bool(common.is_abelian),
                bool(common.is_nilpotent),
                bool(common.is_primitive)
            )

            eigenvalues = common.eigenvalues
            
            db_save(data, eigenvalues)
            print("\n")


