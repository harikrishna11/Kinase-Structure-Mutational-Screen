
def all_motif_groups(motif, offset):

    ret = [("%s_all" % motif, "r%d-%d" % (offset, offset + (len(motif) - 1)))]
    ret.extend([("%s_%d" % (motif, i + 1), "r%d" % (offset + i)) for i in range(len(motif))])
    return ret

def get_rmsd_groups(ape_off, dfg_off, hrd_off, gxgxxg_off, vaik_off, actloop_start, actloop_end, achelix_start, achelix_end, bridge_glut, brige_lys):

    groups = []

    whole_motifs = [("DFG", dfg_off), ("HRD", hrd_off), ("GxGxxG", gxgxxg_off)]
    for (nm, off) in whole_motifs:
        groups.extend(all_motif_groups(nm, off))

    groups.extend([("aC_all", "r%d-%d" % (achelix_start, achelix_end)),
                   ("aC_1h", "r%d-%d" % (achelix_start, (achelix_start + achelix_end) / 2)),
                   ("aC_2h", "r%d-%d" % (((achelix_start + achelix_end) / 2) + 1, achelix_end)),
                   ("aC_start", "r%d" % achelix_start),
                   ("aC_end", "r%d" % achelix_end),
                   ("actloop", "r%d-%d" % (actloop_start, actloop_end)),
                   ("APE_all", "r%d-%d" % (ape_off, ape_off + 2))])

    return groups

def get_dist_groups(ape_off, dfg_off, hrd_off, gxgxxg_off, vaik_off, actloop_start, actloop_end, achelix_start, achelix_end, bridge_glut, brige_lys, dxxxxg_off, rs1_off, rs2_off):

    return [("VAIK_K_NZ", "r%d & a NZ" % (vaik_off + 3)),
            ("bridge_E_CD", "r%d & a CD" % bridge_glut),
            ("DFG_D_CG", "r%d & a CG" % dfg_off),
            ("DFG_F_O", "r%d & a O" % (dfg_off + 1)),
            ("DFG_p2_N", "r%d & a N" % (dfg_off + 4)),
            ("DFG_p1_O", "r%d & a O" % (dfg_off + 3)),
            ("HRD_R_NE", "r%d & a NE" % (hrd_off + 1)),
            #("GxGxxG_5", "r%d" % (gxgxxg_off + 4)), # This is already created for rmsd_groups above, and the only index that uses the dist groups also uses the rmsd groups.
            ("RS1", "r%d" % rs1_off),
            ("RS2", "r%d" % rs2_off),
            ("RS3", "r%d" % (dfg_off + 1)),
            ("RS4", "r%d" % hrd_off),
            ("RSA", "r%d" % dxxxxg_off)]
             

    

    
