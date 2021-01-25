package com.hartwig.hmftools.lilac.read

interface Read {

    fun nucleotide(index: Int): Char

    fun nucleotideIndices(minQual: Int): Collection<Int>

    fun nucleotideIndices(): Collection<Int>

    fun nucleotide(index: Int, minQual: Int): Char

}


