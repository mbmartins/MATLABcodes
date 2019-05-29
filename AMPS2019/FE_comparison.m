%comparison of FE
clear all; close all; clc

SNR =   [60 55  50  45  40  35 30];

mag = load('FE_mag.mat')
ph = load('FE_phase.mat')

plot(SNR, mag.FE_mean, SNR,ph.FE_mean)