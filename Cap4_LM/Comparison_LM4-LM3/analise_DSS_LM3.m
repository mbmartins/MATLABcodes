%
clear all; close all; clc
load("C:\Users\mbrit\Documents\MATLAB\MATLABcodes\Cap4_LM\Comparison_LM4-LM3\comp_LM_DSS")

FE_HLM4_max = max(abs(FE_LM_hist))
FE_DSS_max = max(abs(FE_SS_hist))

TVE_HLM4_max = max(abs(TVE_LM_hist))
TVE_DSS_max = max(abs(TVE_SS_hist))