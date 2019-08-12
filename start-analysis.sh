
MYDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

module load apps/python/2.7.6/gcc-4.4.7

$MYDIR/start-analysis.py "$@"
