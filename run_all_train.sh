#!/bin/bash
if [ "$1" == '' ]; then
  echo "Usage: ./run_all_train.sh iteration"
else
  for i in {1..5};
  do
    ./train $1 model_init.txt "seq_model_0$i.txt" "model_0$i.txt"
  done
fi
