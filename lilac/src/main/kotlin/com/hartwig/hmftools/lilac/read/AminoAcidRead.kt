package com.hartwig.hmftools.lilac.read

interface AminoAcidRead {

    fun aminoAcidIndices(): IntRange

    fun aminoAcid(index: Int, minQual: Int = 0): Char

    fun aminoAcids(startIndex: Int, endIndex: Int, minQual: Int): String {
        val builder = StringBuilder()
        for (i in startIndex..endIndex) {
            val aa = aminoAcid(i, minQual)
            if (aa == '.') {
                return ""
            }

            builder.append(aa)
        }
        return builder.toString()
    }

}


