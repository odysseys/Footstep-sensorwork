#!/bin/sh
#launcher.sh
#navigate to home directory, then to this directory, then execute python script, then back home

cd /
cd home/pi/Documents/Footstep
sleep 10; sudo python readcom_anchor.py
cd /
