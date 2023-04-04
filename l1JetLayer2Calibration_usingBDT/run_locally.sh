#!/bin/bash

#time  python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut --l1MatchGen --MLTarget GenEt            --fracOfDataToUse 0.01 2>&1 | tee cout_ChunkyDonut_GenEt_0.01.txt

#time  python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut --l1MatchGen --MLTarget logGenEt         --fracOfDataToUse 0.01 2>&1 | tee cout_ChunkyDonut_logGenEt_0.01.txt

#time  python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut --l1MatchGen --MLTarget GenEtByL1Et      --fracOfDataToUse 0.01 2>&1 | tee cout_ChunkyDonut_GenEtByL1Et_0.01.txt

#time  python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut --l1MatchGen --MLTarget logGenEtByL1Et   --fracOfDataToUse 0.01 2>&1 | tee cout_ChunkyDonut_logGenEtByL1Et_0.01.txt


time  python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget GenEt            --fracOfDataToUse 0.01 2>&1 | tee cout_PhiRing_GenEt_0.01.txt

time  python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget logGenEt         --fracOfDataToUse 0.01 2>&1 | tee cout_PhiRing_logGenEt_0.01.txt

time  python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget GenEtByL1Et      --fracOfDataToUse 0.01 2>&1 | tee cout_PhiRing_GenEtByL1Et_0.01.txt

time  python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget logGenEtByL1Et   --fracOfDataToUse 0.01 2>&1 | tee cout_PhiRing_logGenEtByL1Et_0.01.txt
