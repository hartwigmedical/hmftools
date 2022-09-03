package com.hartwig.hmftools.cider

import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.genome.region.Strand
import htsjdk.samtools.util.SequenceUtil
import org.apache.logging.log4j.LogManager
import org.eclipse.collections.api.collection.ImmutableCollection
import org.eclipse.collections.api.set.SetIterable
import kotlin.math.ceil
import kotlin.math.floor

fun roundUpToMultiple(v: Int, m: Int) : Int
{
    return (ceil(v / m.toDouble()) * m).toInt()
}

fun roundDownToMultiple(v: Int, m: Int) : Int
{
    return (floor(v / m.toDouble()) * m).toInt()
}

data class AnchorBlosumMatch(
    val anchorStart: Int,
    val anchorEnd: Int,
    val templateAnchorSeq: String,
    val templateGenes: ImmutableCollection<VJAnchorTemplate>,
    val similarityScore: Int)

// create a very simple interface so we can mock the data
interface IAnchorBlosumSearcher
{
    fun searchForAnchor(readString: String) : AnchorBlosumMatch?

    fun searchForAnchor(dnaSeq: String, targetAnchorGeneType: VJGeneType, startOffset: Int = 0,
                        endOffset: Int = dnaSeq.length) : AnchorBlosumMatch?
}

// we want to return a score of whether a sequence
// looks like an anchor
class AnchorBlosumSearcher(
    val ciderGeneDatastore: ICiderGeneDatastore,
    minPartialAnchorAminoAcidLength: Int,
    val allowNegativeSimilarity: Boolean) : IAnchorBlosumSearcher
{
    val minPartialAnchorBaseLength = minPartialAnchorAminoAcidLength * 3

    override fun searchForAnchor(readString: String) : AnchorBlosumMatch?
    {
        var bestMatch: AnchorBlosumMatch? = null

        for (vjGeneType in VJGeneType.values())
        {
            val anchorBlosumMatch: AnchorBlosumMatch? = searchForAnchor(
                readString,
                vjGeneType,
                0,
                readString.length
            )
            if (anchorBlosumMatch != null &&
                (bestMatch == null || anchorBlosumMatch.similarityScore > bestMatch.similarityScore))
            {
                bestMatch = anchorBlosumMatch
            }
        }

        return bestMatch
    }

    override fun searchForAnchor(dnaSeq: String, targetAnchorGeneType: VJGeneType, startOffset: Int, endOffset: Int)
    : AnchorBlosumMatch?
    {
        // sLogger.trace("finding anchor for {}, seq: {}, offset: {}-{}", targetAnchorGeneType, dnaSeq, startOffset, endOffset)

        val templateAnchorSequences : SetIterable<String> = ciderGeneDatastore.getAnchorSequenceSet(targetAnchorGeneType)
        var bestMatch: AnchorBlosumMatch? = null

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
                        if (anchorHomolog != null &&
                            (bestMatch == null || anchorHomolog.similarityScore > bestMatch.similarityScore))
                        {
                            bestMatch = anchorHomolog
                        }
                    }
                }
            }
        }
        else if (targetAnchorGeneType.vj == VJ.V)
        {
            // if we got the V anchor, we have to search in the reverse direction
            // NOTE: this is probably unnecessary
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
                        if (anchorHomolog != null &&
                            (bestMatch == null || anchorHomolog.similarityScore > bestMatch.similarityScore))
                        {
                            bestMatch = anchorHomolog
                        }
                    }
                }
            }
        }

        return bestMatch
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
            val leftTrim: Int = roundUpToMultiple(-anchorStart, 3)
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

            val rightTrim = roundUpToMultiple(anchorEnd - dnaSeq.length, 3)
            anchorEnd -= rightTrim
            trimmedTemplateAnchorSeq = trimmedTemplateAnchorSeq.dropLast(rightTrim)
        }

        val potentialAnchor: String = dnaSeq.substring(anchorStart, anchorEnd)

        assert(trimmedTemplateAnchorSeq.length == potentialAnchor.length)

        if (trimmedTemplateAnchorSeq.length < templateAnchorSeq.length &&
            potentialAnchor.length < minPartialAnchorBaseLength)
        {
            // partial and not enough amino acids
            return null
        }

        if (anchorStart >= 0 && anchorEnd <= dnaSeq.length)
        {
            val score: Int = BlosumSimilarityCalc.calcSimilarityScore(geneType.vj, trimmedTemplateAnchorSeq, potentialAnchor)

            if (score >= 0 || allowNegativeSimilarity)
            {
                val templateGenes: ImmutableCollection<VJAnchorTemplate> = ciderGeneDatastore.getByAnchorSequence(geneType, templateAnchorSeq)
                val anchorBlosumMatch = AnchorBlosumMatch(anchorStart = anchorStart, anchorEnd = anchorEnd,
                    templateAnchorSeq = templateAnchorSeq,
                    templateGenes = templateGenes, similarityScore = score)
                return anchorBlosumMatch
            }
        }
        return null
    }

    companion object
    {
        private val sLogger = LogManager.getLogger(AnchorBlosumSearcher::class.java)
    }
}
