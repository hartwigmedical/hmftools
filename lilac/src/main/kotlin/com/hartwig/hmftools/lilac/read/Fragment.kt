package com.hartwig.hmftools.lilac.read

data class Fragment(val reads: List<Read>) {
    private val aminoAcidIndices = reads
            .flatMap { it.aminoAcidIndices().toList() }
            .distinct()
            .sorted()

    private val nucleotideIndices = reads
            .flatMap { it.nucleotideIndices().toList() }
            .distinct()
            .sorted()


    fun containsAminoAcid(index: Int): Boolean {
        return aminoAcidIndices.contains(index)
    }

    fun aminoAcid(index: Int, minQual: Int = 0): Char {
        for (read in reads) {
            if (read.aminoAcidIndices().contains(index)) {
                return read.aminoAcid(index, minQual)
            }
        }

        throw IllegalArgumentException("Fragment does not contain amino acid at location $index")
    }

    fun nucleotide(index: Int, minQual: Int = 0): Char {
        for (read in reads) {
            if (read.nucleotideIndices().contains(index)) {
                return read.nucleotide(index, minQual)
            }
        }

        throw IllegalArgumentException("Fragment does not contain nucleotide at location $index")
    }

    fun aminoAcidIndices(): List<Int> = aminoAcidIndices

    fun nucleotideIndices(): List<Int> = nucleotideIndices

    fun <T> Iterator<T>.toList(): List<T> {
        val result = mutableListOf<T>()
        while (this.hasNext()) {
            result.add(this.next())
        }
        return result
    }


}