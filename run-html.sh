#!/bin/bash

portNumber="8003"
browserName="wsl-brave"

if [ $# -gt 0 ]; then 
  browserName=$1
fi 

if [ $# -gt 1 ]; then 
  portNumber=$2
fi

if [ $# -gt 2 ]; then 
  echo "ignoring extra arguments to script"
fi

cd payload

if [ $browserName == "firefox" ]; then
  emrun --port $portNumber --browser firefox index.html;
fi 

if [ $browserName == "wsl-edge" ]; then
  emrun --port $portNumber --browser "/mnt/c/Program Files (x86)/Microsoft/Edge/Application/msedge.exe" index.html;
fi 

if [ $browserName == "wsl-brave" ]; then
  emrun --port $portNumber --browser "/mnt/c/Program Files/BraveSoftware/Brave-Browser/Application/brave.exe" index.html;
fi

if [ $browserName == "wsl-chrome" ]; then
  emrun --port $portNumber --browser "/mnt/c/Program Files/Google/Chrome/Application/chrome.exe" index.html;
fi
