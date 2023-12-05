~/.venv\gpu\Scripts\Activate.ps1

$clv_base = 'MYH9_NM_002473.6'
$clv_refseq = 'NM_002473'
$clv_eval = 'MYH9'
$clv_orig = '../clvorig/MYH9.txt'
$clv_variants = ("{0}_clvvariants.csv" -f $clv_base)
$clv_graph = ("{0}_clvgraph.png" -f $clv_base)
../../clvprep.py -o $clv_variants -r $clv_refseq $clv_orig
../../mapvariants.py -o $clv_graph -e $clv_eval -v $clv_variants ("{0}.gb" -f $clv_base)

$dvd_base = $clv_base
$dvd_refseq = $clv_refseq
$dvd_eval = $clv_eval
$dvd_orig = '../dvdorig/MYH9.gvtable.9.csv'
$dvd_variants = ("{0}_dvdvariants.csv" -f $dvd_base)
$dvd_graph = ("{0}_dvdgraph.png" -f $dvd_base)
../../dvdprep.py -o $dvd_variants -r $dvd_refseq $dvd_orig
../../mapvariants.py -o $dvd_graph -e $dvd_eval -v $dvd_variants ("{0}.gb" -f $dvd_base)
