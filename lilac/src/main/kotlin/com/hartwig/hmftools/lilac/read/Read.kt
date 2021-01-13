package com.hartwig.hmftools.lilac.read

interface Read {

    fun nucleotide(index: Int): Char

    fun nucleotideIndices(minQual: Int): Collection<Int>

    fun nucleotideIndices(): Collection<Int>

    fun nucleotide(index: Int, minQual: Int): Char

    fun nucleotides(startIndex: Int, endIndex: Int, minQual: Int): String {
        return sequenceBuilder(startIndex, endIndex) { i:Int -> nucleotide(i, minQual) }
    }

    fun aminoAcidIndices(): Collection<Int>

    fun aminoAcid(index: Int, minQual: Int = 0): Char

    fun aminoAcids(startIndex: Int, endIndex: Int, minQual: Int): String {
        return sequenceBuilder(startIndex, endIndex) { i:Int -> aminoAcid(i, minQual) }
    }

    fun sequenceBuilder(startIndex: Int, endIndex: Int, charFunction: (Int) -> Char): String {
        val builder = StringBuilder()
        for (i in startIndex..endIndex) {
            val aa = charFunction(i)
            if (aa == '.') {
                return ""
            }

            builder.append(aa)
        }
        return builder.toString()
    }

}


