# Reading the sensorvalues of the sensor and send them over serial connection to the computer and save the datas in a file

source: https://learn.adafruit.com/adafruit-10-dof-imu-breakout-lsm303-l3gd20-bmp180/software

All of the libraries in the folder "Libraries" must be installed.

## Linux
Probleme: 
"avrdude: ser_open(): can't open device "/dev/ttyACM0": Permission denied" mit Flatpak; type "ls -l /dev/ttyACM0" into the console and look which group has access to this file. Add your username to this group and restart computer.

## Disable autoreset Arduino on serial connection
https://playground.arduino.cc/Main/DisablingAutoResetOnSerialConnection

