package com.hartwig.hmftools.cider

import com.hartwig.hmftools.common.codon.Codons
import org.apache.logging.log4j.LogManager
import kotlin.math.ceil

fun roundToMultiple(v: Int, m: Int) : Int
{
    return (ceil(v / m.toDouble()) * m).toInt()
}

data class AnchorBlosumMatch(
    val anchorStart: Int,
    val anchorEnd: Int,
    val templateAnchorSeq: String,
    val templateGenes: Collection<VJGene>,
    val similarityScore: Int)

// create a very simple interface so we can mock the data
interface IAnchorBlosumSearcher
{
    fun searchForAnchor(dnaSeq: String, targetAnchorGeneType: VJGeneType, startOffset: Int = 0,
                        endOffset: Int = dnaSeq.length) : AnchorBlosumMatch?
}

// we want to return a score of whether a sequence
// looks like an anchor
class AnchorBlosumSearcher(
    val vjGeneStore: VJGeneStore,
    val minPartialAnchorAminoAcidLength: Int,
    val allowNegativeSimilarity: Boolean) : IAnchorBlosumSearcher
{
    val blosumMapping = BlosumMapping()

    override fun searchForAnchor(dnaSeq: String, targetAnchorGeneType: VJGeneType, startOffset: Int, endOffset: Int)
    : AnchorBlosumMatch?
    {
        sLogger.trace("finding anchor for {}, seq: {}, offset: {}-{}", targetAnchorGeneType, dnaSeq, startOffset, endOffset)

        val templateAnchorSequences : Set<String> = vjGeneStore.getAnchorSequenceSet(targetAnchorGeneType)
        val anchorBlosumMatches = ArrayList<AnchorBlosumMatch>()

        if (targetAnchorGeneType.vj == VJ.J)
        {
            // we want to find the J anchor, by searching forward
            for (i in startOffset until endOffset - 3)
            {
                val aa = Codons.codonToAminoAcid(dnaSeq, i)
                if (aa == CiderUtils.conservedAA(targetAnchorGeneType)) // we might revisit this case later
                {
                    // this is a potential match
                    for (templateAnchorSeq in templateAnchorSequences)
                    {
                        val start: Int = i
                        val end: Int = i + templateAnchorSeq.length
                        // make sure the anchor is long enough, for now we don't allow short anchors
                        val anchorHomolog = tryMatchWithBlosum(targetAnchorGeneType, dnaSeq, start, end, templateAnchorSeq)
                        if (anchorHomolog != null && (allowNegativeSimilarity || anchorHomolog.similarityScore >= 0))
                            anchorBlosumMatches.add(anchorHomolog)
                    }
                }
            }
        }
        else if (targetAnchorGeneType.vj == VJ.V)
        {
            // if we got the V anchor, we have to search in the reverse direction
            for (i in (endOffset - 1) downTo  startOffset + 3)
            {
                val aa = Codons.codonToAminoAcid(dnaSeq, i - 3)
                if (aa == CiderUtils.conservedAA(targetAnchorGeneType))
                {
                    // V anchor ends with C
                    for (templateAnchorSeq in templateAnchorSequences)
                    {
                        val start: Int = i - templateAnchorSeq.length
                        val end: Int = i
                        val anchorHomolog = tryMatchWithBlosum(targetAnchorGeneType, dnaSeq, start, end, templateAnchorSeq)
                        if (anchorHomolog != null && (allowNegativeSimilarity || anchorHomolog.similarityScore >= 0))
                            anchorBlosumMatches.add(anchorHomolog)
                    }
                }
            }
        }

        if (anchorBlosumMatches.isEmpty())
            return null

        // sort them by most similar to least
        anchorBlosumMatches.sortByDescending({ o -> o.similarityScore })

        return anchorBlosumMatches.maxByOrNull({ match -> match.similarityScore })
    }

    private fun tryMatchWithBlosum(
        geneType: VJGeneType,
        dnaSeq: String,
        inputAnchorStart: Int,
        inputAnchorEnd: Int,
        templateAnchorSeq: String) : AnchorBlosumMatch?
    {
        var anchorStart = inputAnchorStart
        var anchorEnd = inputAnchorEnd
        var trimmedTemplateAnchorSeq = templateAnchorSeq

        // we might want to deal with partial sequence matches, by trimming the
        // template anchor to be same size as the input anchor
        if (anchorStart < 0)
        {
            if (geneType.vj == VJ.J)
            {
                // we don't want to allow removing the conserved F or W AA
                return null
            }

            // since we are searching backwards in the sequence, we just have to deal
            // with partial anchor match, and we want to make sure the full length
            // is multiple of 3
            val leftTrim: Int = roundToMultiple(-anchorStart, 3)
            anchorStart += leftTrim
            trimmedTemplateAnchorSeq = trimmedTemplateAnchorSeq.substring(leftTrim)
        }

        if (anchorEnd >= dnaSeq.length)
        {
            if (geneType.vj == VJ.V)
            {
                // we don't want to allow removing the conserved C
                return null
            }

            val rightTrim = roundToMultiple(anchorEnd - dnaSeq.length, 3)
            anchorEnd -= rightTrim
            trimmedTemplateAnchorSeq = trimmedTemplateAnchorSeq.dropLast(rightTrim)
        }

        val templateAnchorAA: String = Codons.aminoAcidFromBases(trimmedTemplateAnchorSeq)

        val potentialAnchor: String = dnaSeq.substring(anchorStart, anchorEnd)

        assert(trimmedTemplateAnchorSeq.length == potentialAnchor.length)

        val potentialAnchorAA: String = Codons.aminoAcidFromBases(potentialAnchor)

        if (trimmedTemplateAnchorSeq.length < templateAnchorSeq.length &&
            potentialAnchorAA.length < minPartialAnchorAminoAcidLength)
        {
            // partial and not enough amino acids
            return null
        }

        if (anchorStart >= 0 && anchorEnd <= dnaSeq.length)
        {
            if (potentialAnchorAA.contains(Codons.UNKNOWN))
                return null

            val score: Int = calcSimilarityScore(templateAnchorAA, potentialAnchorAA)

            //if (score >= 0)
            //{
                /*sLogger.debug(
                    "anchor found seq: {}, anchor type: {}, anchor seq: {}",
                    dnaSeq, geneType, potentialAnchor
                )*/

                val templateGenes: Collection<VJGene> = vjGeneStore.getByAnchorSequence(geneType, templateAnchorSeq)
                val anchorBlosumMatch = AnchorBlosumMatch(anchorStart = anchorStart, anchorEnd = anchorEnd,
                    templateAnchorSeq = templateAnchorSeq,
                    templateGenes = templateGenes, similarityScore = score)
                return anchorBlosumMatch
            //}
        }
        return null
    }

    // we do with amino acid sequence
    // diff of 0 means exact match
    fun calcBlosumDistance(refAnchorAA: String, seqAA: String) : Int
    {
        if (refAnchorAA.length != seqAA.length)
            return Int.MAX_VALUE

        // calculate the blosum score
        return blosumMapping.calcSequenceSum(refAnchorAA) - blosumMapping.calcSequenceSum(refAnchorAA, seqAA)
    }

    fun calcSimilarityScore(refAnchorAA: String, seqAA: String) : Int
    {
        return CiderConstants.MAX_BLOSUM_DIFF_PER_AA * seqAA.length -
                CiderConstants.BLOSUM_SIMILARITY_SCORE_CONSTANT -
                calcBlosumDistance(refAnchorAA, seqAA)
    }

    companion object
    {
        private val sLogger = LogManager.getLogger(AnchorBlosumSearcher::class.java)
    }
}
