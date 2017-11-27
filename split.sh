#!/bin/bash

samtools view $1 | python /src/scour.py > $2 2> $3
