CCX="g++"

set -e

mkdir -p build/

set -x

$CCX -O3 -march=native src/exp/rng_throughput.cpp -o build/rng_throughput
$CCX -O3 -march=native src/exp/hqe_scheme_prep.cpp -o build/hqe_scheme_prep

