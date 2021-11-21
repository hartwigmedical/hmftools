package com.hartwig.hmftools.teal.util

import com.hartwig.hmftools.common.aligner.SequenceAligner
import com.hartwig.hmftools.teal.TeloConstants
import org.apache.commons.lang3.StringUtils
import org.apache.logging.log4j.LogManager

// helper class to match a sequence to telomere
// and give a similarity rate
object TelomereMatcher
{
    private val LOGGER = LogManager.getLogger(javaClass)

    data class TelomereMatch(val start: Int, val numCompleteMatch: Int, val sequence: String)

    fun findGTelomereSegment(seq: String, matchThreshold: Double): TelomereMatch?
    {
        return findTelomereSegmentImpl(seq, matchThreshold, gRich = true)
    }

    fun findCTelomereSegment(seq: String, matchThreshold: Double): TelomereMatch?
    {
        return findTelomereSegmentImpl(seq, matchThreshold, gRich = false)
    }

    // calculate how much this sequence resembles telomere
    private fun findTelomereSegmentImpl(seq: String, matchThreshold: Double, gRich: Boolean): TelomereMatch?
    {
        // we want to allow some insertion / deletion plus a full hexamer
        val telomereTemplateLength = (seq.length * 1.2).toInt() + 6

        // we provide a template that is 1.5x as long as the original sequence
        val alignOps = SequenceAligner.alignSubsequence(seq, generateTelomereTemplate(telomereTemplateLength, gRich))

        // we work out longest stretch of telomere which satisfy our matching threshold

        val scoreSeq = alignOps.map { op -> if (op == SequenceAligner.AlignOp.MATCH) 1.0 else 0.0 }.toDoubleArray()

        val longestMatchRange: IntRange = LongestSegment.longestSegmentAverage(scoreSeq, matchThreshold)

        if (longestMatchRange.isEmpty())
            return null

        var matchStart = longestMatchRange.first
        var matchEnd = longestMatchRange.last + 1

        // we want to trim the ends, cause the longest stretch could well include some mismatches at either ends
        while (matchStart > matchEnd && alignOps[matchStart] != SequenceAligner.AlignOp.MATCH)
        {
            ++matchStart
        }

        while (matchStart > matchEnd && alignOps[matchEnd - 1] != SequenceAligner.AlignOp.MATCH)
        {
            --matchEnd
        }

        // now we go backwards and work out the match sequence
        val seqBuilder = StringBuilder()
        val seqItr = seq.iterator()
        var numMatch = 0
        val matchSegmentLength = matchEnd - matchStart
        for (i in alignOps.indices)
        {
            val alignOp = alignOps[i]
            when (alignOp)
            {
                SequenceAligner.AlignOp.MATCH ->
                {
                    if (i >= matchStart && i < matchEnd)
                    {
                        ++numMatch
                    }
                }
                SequenceAligner.AlignOp.DELETION ->
                {
                    continue
                }
                else -> {}
            }
            val b = seqItr.next()
            if (i in matchStart until matchEnd)
            {
                seqBuilder.append(b)
            }
        }
        val matchRatio = numMatch.toDouble() / matchSegmentLength
        val matchSeq = seqBuilder.toString()

        LOGGER.trace("seq({}) matchSegmentLength({}) numMatch({}) matchSegment({}) ratio({})", seq, matchSegmentLength, numMatch, matchSeq, matchRatio)
        return TelomereMatch(matchStart, numMatch, matchSeq)
    }

    private fun generateTelomereTemplate(length: Int, gRich: Boolean): String
    {
        val canonicalHexamer = if (gRich) TeloConstants.CANONICAL_TELOMERE_SEQ else TeloConstants.CANONICAL_TELOMERE_SEQ_REV
        val b = StringBuilder()
        val numFullHexamers = length / canonicalHexamer.length
        b.append(StringUtils.repeat(canonicalHexamer, numFullHexamers))
        val residualLength = length % canonicalHexamer.length
        if (residualLength > 0)
        {
            b.append(canonicalHexamer, 0, residualLength)
        }
        return b.toString()
    }

    fun matchesGTelomere(seq: String, telomereMatchThreshold: Double, minTelomereMatchLength: Int): Boolean
    {
        if (seq.length < minTelomereMatchLength)
            return false

        if (!seq.contains("GGG"))
        {
            return false
        }

        val telomereMatch = findGTelomereSegment(seq, telomereMatchThreshold)
        return telomereMatch != null && telomereMatch.numCompleteMatch >= minTelomereMatchLength
    }

    fun matchesCTelomere(seq: String, telomereMatchThreshold: Double, minTelomereMatchLength: Int): Boolean
    {
        if (seq.length < minTelomereMatchLength)
            return false

        if (!seq.contains("CCC"))
        {
            return false
        }

        val telomereMatch = findCTelomereSegment(seq, telomereMatchThreshold)
        return telomereMatch != null && telomereMatch.numCompleteMatch >= minTelomereMatchLength
    }
}
