
import _mysql
import os
import os.path

db = _mysql.connect(host = "acbbdb1.picr.man.ac.uk", user = "nsmd", passwd = "nsmdP123", db = "nsmd")
data_root = "/data/snc/nstephenson/nsmd/results"

dirs = os.listdir(data_root)

def is_numerical_dir(name):

    fullname = os.path.join(data_root, name)
    if not os.path.isdir(fullname):
        return False
    try:
        n = int(name)
        return True
    except ValueError:
        return False

dirs = filter(is_numerical_dir, dirs)

for d in dirs:
    
    db.query("select gene, mutation from sequences inner join kinases on kinaseid = kinases.id inner join genes on gene_id = genes.id where sequences.id = %s" % d)
    rows = db.store_result().fetch_row(maxrows=0)

    if len(rows) == 0:
        print >>sys.stderr, "Ignore directory", d, "not in DB"
        continue
    elif len(rows) > 1:
        raise Exception("Duplicate entries in DB for %s" % d)

    protdir = os.path.join(data_root, "by-protein", rows[0][0])
    try:
        os.makedirs(protdir)
    except OSError as e:
        if e.errno != 17: # Already exists
            raise e
    
    linkname = os.path.join(protdir, rows[0][1].replace(" ", "_"))

    try:
        os.symlink(os.path.join(data_root, d), linkname)
    except OSError as e:
        if e.errno != 17:
            raise e

