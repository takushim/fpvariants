~/.venv\gpu\Scripts\Activate.ps1

$basename = 'MYO3A_NM_017433.5'
$eval = 'MYO3A'
$output_type = 'png'

$clv_orig = '../clvorig/MYO3A.txt'
$clv_variants = ("{0}_clvvariants.csv" -f $basename)
$clv_graph = ("{0}_clvgraph.{1}" -f $basename, $output_type)
../../clvprep.py -o $clv_variants -e $eval $clv_orig
../../mapvariants.py -o $clv_graph -e $eval -v $clv_variants ("{0}.gb" -f $basename)

$dvd_orig = '../dvdorig/MYO3A.gvtable.9.csv'
$dvd_variants = ("{0}_dvdvariants.csv" -f $basename)
$dvd_graph = ("{0}_dvdgraph.{1}" -f $basename, $output_type)
#../../dvdprep.py -o $dvd_variants -e $eval $dvd_orig
#../../mapvariants.py -o $dvd_graph -e $eval -v $dvd_variants ("{0}.gb" -f $basename)
