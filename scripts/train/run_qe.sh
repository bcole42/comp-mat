#!/usr/bin/env bash
# Run pw.x on every training/validation input. Usage:  NP=4 bash 02_run_qe.sh
NP=${NP:-4}

cd qe || exit 1
for inp in si_*.in; do
  out="${inp%.in}.out"
  if [ -f "$out" ] && grep -q "JOB DONE" "$out"; then echo "skip $inp"; continue; fi
  echo "run  $inp"
  mpirun -np "$NP" pw.x -in "$inp" > "$out"
done
cd ..
echo "Done."
