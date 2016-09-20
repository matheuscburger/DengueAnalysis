#!/bin/bash

src/preprocess_Illumina.sh &> log/preprocess_Illumina.txt
src/preprocess_Affymetrix.sh &> log/preprocess_Affymetrix.txt
