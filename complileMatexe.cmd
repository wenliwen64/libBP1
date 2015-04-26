#!/bin/bash -f
qrsh -l i,h_rt=4:00:00
 module load matlab
mcc -m runteleBPcontHoff.m
