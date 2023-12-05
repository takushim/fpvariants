~/.venv\gpu\Scripts\Activate.ps1

$basename = 'MYO3A_NM_017433.5'

$clv_orig = '../clvorig/MYO3A.txt'
$clv_variants = ("{0}_clvvariants.csv" -f $basename)
$clv_graph = ("{0}_clvgraph.png" -f $basename)
../../clvprep.py -o $clv_variants $clv_orig
../../mapvariants.py -o $clv_graph -e MYO3A -v $clv_variants ("{0}.gb" -f $basename)

$dvd_orig = '../dvdorig/MYO3A.gvtable.9.csv'
$dvd_variants = ("{0}_dvdvariants.csv" -f $basename)
$dvd_graph = ("{0}_dvdgraph.png" -f $basename)
../../dvdprep.py -o $dvd_variants $dvd_orig
../../mapvariants.py -o $dvd_graph -e MYO3A -v $dvd_variants ("{0}.gb" -f $basename)
