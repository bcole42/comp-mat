#!/usr/bin/env bash
# Run pw.x on every input in conv_ecut/ and conv_k/.
# Usage:  NP=4 bash 02_run_all.sh          (NP = number of MPI processes)

NP=${NP:-4}

for dir in conv_ecut conv_k; do
  cd "$dir" || exit 1
  for inp in *.in; do
    out="${inp%.in}.out"
    if [ -f "$out" ] && grep -q "JOB DONE" "$out"; then
      echo "skip  $dir/$inp (already done)"
      continue
    fi
    echo "run   $dir/$inp"
    mpirun -np "$NP" pw.x -in "$inp" > "$out"
  done
  cd ..
done

echo "Done."
