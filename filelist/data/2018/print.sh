#!/bin/bash
printf '%s\n' "/hadoop/cms/store/user/legianni/skimNanoVVH-Hadronic__v2/JetHT/Run2018A-UL2018_MiniAODv2_NanoAODv9-v2/"*"/0000"/*.root > data_2018A.txt
printf '%s\n' "/hadoop/cms/store/user/legianni/skimNanoVVH-Hadronic__v2/JetHT/Run2018B-UL2018_MiniAODv2_NanoAODv9-v1/"*"/0000"/*.root > data_2018B.txt
printf '%s\n' "/hadoop/cms/store/user/legianni/skimNanoVVH-Hadronic__v2/JetHT/Run2018C-UL2018_MiniAODv2_NanoAODv9-v1/"*"/0000"/*.root > data_2018C.txt
printf '%s\n' "/hadoop/cms/store/user/legianni/skimNanoVVH-Hadronic__v2/JetHT/Run2018D-UL2018_MiniAODv2_JMENanoAODv9-v1/"*"/0000"/*.root > data_2018D.txt