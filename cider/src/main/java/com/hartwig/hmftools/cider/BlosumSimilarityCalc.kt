package com.hartwig.hmftools.cider

import com.hartwig.hmftools.common.codon.Codons
import org.apache.logging.log4j.LogManager

object BlosumSimilarityCalc
{
    private val sLogger = LogManager.getLogger(javaClass)
    val blosumMapping = BlosumMapping()

    // we do with amino acid sequence
    // diff of 0 means exact match
    private fun calcBlosumDistance(refAnchorAA: String, seqAA: String) : Int
    {
        if (refAnchorAA.length != seqAA.length)
            return Int.MAX_VALUE

        try
        {
            // calculate the blosum score
            return blosumMapping.calcSequenceSum(refAnchorAA) - blosumMapping.calcSequenceSum(refAnchorAA, seqAA)
        }
        catch (e: IllegalArgumentException)
        {
            throw IllegalArgumentException("cannot calc blosum for seq: ${refAnchorAA} and ${seqAA}")
        }
    }

    private fun calcAminoAcidSimilarityScore(refAA: String, seqAA: String) : Int
    {
        val similarity = CiderConstants.MAX_BLOSUM_DIFF_PER_AA * refAA.length -
                CiderConstants.BLOSUM_SIMILARITY_SCORE_CONSTANT -
                calcBlosumDistance(refAA, seqAA)

        //sLogger.trace("{} vs {} similiarity={}", refAA, seqAA, similarity)

        return similarity
    }

    @JvmStatic
    fun calcSimilarityScore(vj: VJ, templateDnaSeq: String, dnaSeq: String) : Int
    {
        var templateDnaSeqFixed = templateDnaSeq
        var dnaSeqFixed = dnaSeq

        // trim them to same length
        if (templateDnaSeq.length != dnaSeq.length)
        {
            // for V we must align to right side, and set length to multiple of 3
            val trimToLength = roundDownToMultiple(Math.min(templateDnaSeq.length, dnaSeq.length), 3)

            if (vj == VJ.V)
            {
                // align right
                templateDnaSeqFixed = templateDnaSeqFixed.takeLast(trimToLength)
                dnaSeqFixed = dnaSeqFixed.takeLast(trimToLength)
            }
            else
            {
                // align left
                templateDnaSeqFixed = templateDnaSeqFixed.take(trimToLength)
                dnaSeqFixed = dnaSeqFixed.take(trimToLength)
            }
        }

        // try to fix up for cases where sequence contains N
        if (templateDnaSeqFixed.contains('N') || dnaSeqFixed.contains('N'))
        {
            val templateDnaSeqBuilder = StringBuilder(templateDnaSeqFixed)
            val dnaSeqBuilder = StringBuilder(dnaSeqFixed)

            // change any unknown base N to the other one
            for (i in templateDnaSeqFixed.indices)
            {
                if (templateDnaSeqFixed[i] == 'N')
                {
                    templateDnaSeqBuilder[i] = dnaSeqFixed[i]
                }
                if (dnaSeqFixed[i] == 'N')
                {
                    dnaSeqBuilder[i] = templateDnaSeqFixed[i]
                }
            }

            templateDnaSeqFixed = templateDnaSeqBuilder.toString()
            dnaSeqFixed = templateDnaSeqBuilder.toString()

            if (templateDnaSeqFixed.contains('N') || dnaSeqFixed.contains('N'))
            {
                throw IllegalArgumentException("calcSimilarityScore: both (${templateDnaSeq} and ${dnaSeq}) contain N at same location")
            }
        }

        val templateAaSeq = Codons.aminoAcidFromBases(templateDnaSeqFixed)
        val aaSeq = Codons.aminoAcidFromBases(dnaSeqFixed)

        return calcAminoAcidSimilarityScore(templateAaSeq, aaSeq)
    }
}