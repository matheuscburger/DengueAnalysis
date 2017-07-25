#!/bin/bash

echo "Preprocessing Illumina ..."
src/preprocess_Illumina.sh &> log/preprocess/preprocess_Illumina.txt || { echo "Unable to preprocess Illumina data"; exit 1; }
echo "Preprocessing Affymetrix ..."
src/preprocess_Affymetrix.sh &> log/preprocess/preprocess_Affymetrix.txt || { echo "Unable to preprocess Affymetrix data"; exit 1; }
echo "Finalizing ... "
src/final_preprocess.sh &> log/preprocess/final_preprocess.txt || { echo "Unable to preprocess data"; exit 1; }
echo "Preprocess done."
