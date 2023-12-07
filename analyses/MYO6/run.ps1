~/.venv\gpu\Scripts\Activate.ps1


$clv_base = 'MYO6_NM_004999.4'
$clv_refseq = 'NM_004999'
$clv_eval = 'MYO6'
$clv_orig = '../clvorig/MYO6.txt'
$clv_variants = ("{0}_clvvariants.csv" -f $clv_base)
$clv_graph = ("{0}_clvgraph.pdf" -f $clv_base)
../../clvprep.py -o $clv_variants -r $clv_refseq -e $clv_eval $clv_orig
../../mapvariants.py -o $clv_graph -e $clv_eval -v $clv_variants -s 0 4500 ("{0}.gb" -f $clv_base)

#$dvd_base = 'MYO6_NM_001368865.1'
#$dvd_refseq = 'NM_001368865'
#$dvd_eval = $clv_eval
#$dvd_orig = '../dvdorig/MYO6.gvtable.9.csv'
#$dvd_variants = ("{0}_dvdvariants.csv" -f $dvd_base)
#$dvd_graph = ("{0}_dvdgraph.png" -f $dvd_base)
#../../dvdprep.py -o $dvd_variants -r $dvd_refseq -e $clv_eval $dvd_orig
#../../mapvariants.py -o $dvd_graph -e $dvd_eval -v $dvd_variants -s 0 4500 ("{0}.gb" -f $dvd_base)
