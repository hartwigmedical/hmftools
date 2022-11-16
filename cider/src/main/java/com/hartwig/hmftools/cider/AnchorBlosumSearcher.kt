package com.hartwig.hmftools.cider

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
    enum class Mode
    {
        ALLOW_NEG_SIMILARITY,
        DISALLOW_NEG_SIMILARITY
    }

    fun searchForAnchor(readString: String, mode: Mode) : AnchorBlosumMatch?

    fun searchForAnchor(sequence: String,
                        targetAnchorGeneTypes: Collection<VJGeneType>,
                        mode: Mode,
                        startOffset: Int = 0, // search start offset in sequence
                        endOffset: Int = sequence.length // search end offset in sequence
    ) : AnchorBlosumMatch?
}

// we want to return a score of whether a sequence
// looks like an anchor
class AnchorBlosumSearcher(
    val ciderGeneDatastore: ICiderGeneDatastore,
    minPartialAnchorAminoAcidLength: Int) : IAnchorBlosumSearcher
{
    val minPartialAnchorBaseLength = minPartialAnchorAminoAcidLength * 3

    override fun searchForAnchor(readString: String, mode: IAnchorBlosumSearcher.Mode) : AnchorBlosumMatch?
    {
        return searchForAnchor(
            readString,
            VJGeneType.values().toList(),
            mode,
            0,
            readString.length
        )
    }

    override fun searchForAnchor(sequence: String,
                                 targetAnchorGeneTypes: Collection<VJGeneType>,
                                 mode: IAnchorBlosumSearcher.Mode,
                                 startOffset: Int,
                                 endOffset: Int)
    : AnchorBlosumMatch?
    {
        // sLogger.trace("finding anchor for {}, seq: {}, offset: {}-{}", targetAnchorGeneType, dnaSeq, startOffset, endOffset)

        var bestMatch: AnchorBlosumMatch? = null

        for (targetAnchorGeneType in targetAnchorGeneTypes)
        {
            val templateAnchorSequences : SetIterable<String> = ciderGeneDatastore.getAnchorSequenceSet(targetAnchorGeneType)

            // We match each template anchor against the input DNA
            for (i in startOffset  until endOffset)
            {
                //val aa = Codons.codonToAminoAcid(dnaSeq, i)

                //if (aa == CiderUtils.conservedAA(targetAnchorGeneType)) // we might revisit this case later
                //{
                    // this is a potential match
                    for (templateAnchorSeq in templateAnchorSequences)
                    {
                        val anchorPos: Int

                        if (targetAnchorGeneType.vj == VJ.V)
                        {
                            // to allow for partial sequences, for V type anchor, i is the last base of
                            // the anchor
                            // ------------------------------------------- input DNA sequence
                            //         ++++++++++    template anchor
                            //                  |
                            //                  i
                            anchorPos = i - templateAnchorSeq.length + 1
                        }
                        else
                        {
                            // to allow for partial sequences, for J type anchor, i is the first base of
                            // the anchor
                            // ------------------------------------------- input DNA sequence
                            //         ++++++++++    template anchor
                            //         |
                            //         i
                            anchorPos = i
                        }

                        // make sure the anchor is long enough, for now we don't allow short anchors
                        val anchorHomolog = tryMatchWithBlosum(targetAnchorGeneType, sequence, anchorPos,
                            templateAnchorSeq, mode)

                        if (anchorHomolog != null &&
                            (bestMatch == null || anchorHomolog.similarityScore > bestMatch.similarityScore))
                        {
                            bestMatch = anchorHomolog
                        }
                    }
                //}
            }
        }

        return bestMatch
    }

    private fun tryMatchWithBlosum(
        geneType: VJGeneType,
        dnaSeq: String,
        inputAnchorStart: Int,
        templateAnchorSeq: String,
        mode: IAnchorBlosumSearcher.Mode) : AnchorBlosumMatch?
    {
        var anchorStart = inputAnchorStart
        var anchorEnd = inputAnchorStart + templateAnchorSeq.length
        var trimmedTemplateAnchorSeq = templateAnchorSeq

        // we want to deal with partial sequence matches, by trimming the
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

            if (score >= 0 || mode == IAnchorBlosumSearcher.Mode.ALLOW_NEG_SIMILARITY)
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
