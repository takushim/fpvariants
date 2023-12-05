~/.venv\gpu\Scripts\Activate.ps1

$variant_file = '../clinvar/MYO15A.txt'
$basename = 'test'

../clvprep.py -o ("{0}_clvariants.csv" -f $basename) $variant_file

#. ../dvdmap.py -m 290 MYO7A_NM_000260.4.gb