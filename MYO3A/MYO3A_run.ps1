~/.venv\gpu\Scripts\Activate.ps1

$variant_file = '../clinvar/MYO3A.txt'
$basename = 'MYO3A_NM_017433.5'

../clvprep.py -o ("{0}_clvariants.csv" -f $basename) $variant_file

