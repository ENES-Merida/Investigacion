#!/bin/bash
#nvc -fast -acc=gpu -gpu=cc75 -Minfo=all *.c -o test
nvc -acc=gpu -gpu=cc75 -Minfo=all $1 -o test
