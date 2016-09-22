#!/bin/bash

echo "Preprocessing Illumina ..."
src/preprocess_Illumina.sh &> log/preprocess_Illumina.txt
echo "Preprocessing Affymetrix ..."
src/preprocess_Affymetrix.sh &> log/preprocess_Affymetrix.txt
echo "Finalizing ... "
src/final_preprocess.sh &> log/final_preprocess.txt
echo "Preprocess done."
