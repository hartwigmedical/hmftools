package com.hartwig.hmftools.lilac.read

data class Fragment(val reads: List<AminoAcidRead>) {
    private val aminoAcidIndices = reads
            .flatMap { it.aminoAcidIndices().toList() }
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

    fun aminoAcidIndices(): List<Int> = aminoAcidIndices


    fun <T> Iterator<T>.toList(): List<T> {
        val result = mutableListOf<T>()
        while (this.hasNext()) {
            result.add(this.next())
        }
        return result
    }


}