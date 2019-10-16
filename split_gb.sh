#!/bin/bash

awk -v n=1 '/^\/\//{close("out"n);n++;next} {print > "out"n}'

