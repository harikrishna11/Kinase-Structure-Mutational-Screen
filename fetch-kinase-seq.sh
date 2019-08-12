#!/bin/bash

echo "select MID(canonical_sequence, offset_canonical, length) from genes join kinases on gene_id = genes.id where gene = \"$1\"" | mysql -h acbbdb1 -u nsmd --password=nsmdP123 nsmd | tail -n +2
