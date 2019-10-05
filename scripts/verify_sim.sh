#!/bin/bash
intersect "00.bed" $1/present.bed -v -f 0.8 > 00_as_00.bed
intersect "00.bed" $1/heterozygous.bed -f 0.8 > 10_as_00.bed
intersect "00.bed" $1/homozygous.bed -f 0.8 > 11_as_00.bed

intersect "10.bed" $1/present.bed -v -f 0.8 > 00_as_10.bed
intersect "10.bed" $1/heterozygous.bed -f 0.8 > 10_as_10.bed
intersect "10.bed" $1/homozygous.bed -f 0.8 > 11_as_10.bed

intersect "11.bed" $1/present.bed -v -f 0.8 > 00_as_11.bed
intersect "11.bed" $1/heterozygous.bed -f 0.8 > 10_as_11.bed
intersect "11.bed" $1/homozygous.bed -f 0.8 > 11_as_11.bed
