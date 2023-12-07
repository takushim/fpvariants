~/.venv\gpu\Scripts\Activate.ps1

$basename = 'MYO7A_NM_000260.4'
$eval = 'MYO7A'
$refseq = 'NM_000260'

$clv_orig = '../clvorig/MYO7A.txt'
$clv_variants = ("{0}_clvvariants.csv" -f $basename)
$clv_graph = ("{0}_clvgraph.pdf" -f $basename)
../../clvprep.py -o $clv_variants -r $refseq -e $eval $clv_orig
../../mapvariants.py -o $clv_graph -e $eval -v $clv_variants ("{0}.gb" -f $basename)

#$dvd_orig = '../dvdorig/MYO7A.gvtable.9.csv'
#$dvd_variants = ("{0}_dvdvariants.csv" -f $basename)
#$dvd_graph = ("{0}_dvdgraph.png" -f $basename)
#../../dvdprep.py -o $dvd_variants -r $refseq -e $eval $dvd_orig
#../../mapvariants.py -o $dvd_graph -e $eval -v $dvd_variants ("{0}.gb" -f $basename)
