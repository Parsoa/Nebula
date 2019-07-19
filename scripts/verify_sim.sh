#!/bin/bash
p="00"
intersect "00.bed" $1/absent.bed > 00_as_00.bed
intersect "00.bed" $1/heterozygous.bed > 10_as_00.bed
intersect "00.bed" $1/homozygous.bed > 11_as_00.bed
p="10"
intersect "10.bed" $1/absent.bed > 00_as_10.bed
intersect "10.bed" $1/heterozygous.bed > 10_as_10.bed
intersect "10.bed" $1/homozygous.bed > 11_as_10.bed
p="11"
intersect "11.bed" $1/absent.bed > 00_as_11.bed
intersect "11.bed" $1/heterozygous.bed > 10_as_11.bed
intersect "11.bed" $1/homozygous.bed > 11_as_11.bed
