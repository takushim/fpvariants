~/.venv\gpu\Scripts\Activate.ps1

$variant_file = '../clinvar/MYO7A.txt'
$basename = 'MYO7A_NM_000260.4'

../clvprep.py -o ("{0}_clinvar.csv" -f $basename) $variant_file

#. ../dvdmap.py -m 290 MYO7A_NM_000260.4.gb