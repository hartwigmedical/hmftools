package htsjdk.samtools

object SAMUtilsExtra {

    fun SAMRecord.getAlignmentBlocksWithSoftClips(): List<AlignmentBlock> {
        return getAlignmentBlocksWithSoftClips(cigar, alignmentStart, "read cigar")
    }

    fun getAlignmentBlocksWithSoftClips(cigar: Cigar, alignmentStart: Int, cigarTypeName: String): List<AlignmentBlock> {
        val alignmentBlocks= mutableListOf<AlignmentBlock>()
        var readBase = 1
        var refBase = alignmentStart
        for (i in 0 until cigar.cigarElements.size) {
            val current = cigar.cigarElements[i]!!

            when (current.operator) {
                CigarOperator.H -> {
                }
                CigarOperator.P -> {
                }
                CigarOperator.S-> {
                    // INCLUDE FRONT SOFT CLIP
                    if (i == 0 && i < cigar.cigarElements.size - 1) {
                        val next = cigar.cigarElements[i + 1].operator!!
                        if (next  == CigarOperator.M) {
                            alignmentBlocks.add(AlignmentBlock(readBase, refBase - current.length, current.length))
                        }
                    }

                    // INCLUDE END SOFT CLIP
                    if (i == cigar.cigarElements.size - 1 && i > 0) {
                        val prev = cigar.cigarElements[i - 1].operator!!
                        if (prev == CigarOperator.M) {
                            alignmentBlocks.add(AlignmentBlock(readBase, refBase, current.length))
                        }
                    }

                    readBase += current.length
                }
                CigarOperator.N -> refBase += current.length
                CigarOperator.D -> refBase += current.length
                CigarOperator.I -> readBase += current.length
                CigarOperator.M, CigarOperator.EQ, CigarOperator.X -> {
                    val length = current.length
                    alignmentBlocks.add(AlignmentBlock(readBase, refBase, length))
                    readBase += length
                    refBase += length
                }
                else -> throw IllegalStateException("Case statement didn't deal with " + cigarTypeName + " op: " + current.operator + "in CIGAR: " + cigar)
            }

        }


        return alignmentBlocks
    }

}