~/.venv\gpu\Scripts\Activate.ps1

$folders = @("MYH9", "MYH14", "MYO3A", "MYO6", "MYO7A", "MYO15A")

foreach ($folder in $folders) {
    cd $folder
    ./run.ps1
    cd ..
}