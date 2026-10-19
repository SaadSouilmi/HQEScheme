CCX="g++"

mkdir -p build/
set -x

$CCX -O3 -march=native src/exp/rng_throughput.cpp -o build/rng_throughput

set +x
