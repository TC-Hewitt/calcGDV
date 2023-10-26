#!/usr/bin/awk -f

# Usage: ./transpose.awk in.txt > out.txt

awk
    BEGIN { FS = OFS = "\t" }
    !/^#/ {
        for (i = 1; i <= NF; i++) {
            a[i,NR] = $i
            if (max_nf < NF) max_nf = NF
        }
    }
    END {
        for (i = 1; i <= max_nf; i++) {
            for (j = 1; j <= NR; j++) {
                printf "%s%s", a[i,j], (j == NR ? "\n" : OFS)
            }
        }
    }