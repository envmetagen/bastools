#!/bin/bash

awk -v n=1 '/^\/\//{close("outTemp"n);n++;next} {print > "outTemp"n}'

