package com.hartwig.hmftools.teal.util

import com.hartwig.hmftools.common.sequence.*
import com.hartwig.hmftools.teal.TealConstants
import org.apache.logging.log4j.LogManager

// helper class to match a sequence to telomere
// and give a similarity rate
object TelomereMatcher
{
    private val LOGGER = LogManager.getLogger(javaClass)

    data class TelomereMatch(val matchStart: Int, val matchEnd: Int, val numMatchedBases: Int, val matchedSequence: String)

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

        val longestMatchRange = LongestSegment.longestSegmentAverage(scoreSeq, matchThreshold)

        if (longestMatchRange == null)
        {
            return null
        }

        var alignOpStart = longestMatchRange.minimum
        var alignOpEnd = longestMatchRange.maximum + 1

        // we want to trim the ends, cause the longest stretch could well include some mismatches at either ends
        while (alignOpStart < alignOpEnd && alignOps[alignOpStart] != SequenceAligner.AlignOp.MATCH)
        {
            ++alignOpStart
        }

        while (alignOpStart < alignOpEnd && alignOps[alignOpEnd - 1] != SequenceAligner.AlignOp.MATCH)
        {
            --alignOpEnd
        }

        if (alignOpStart == alignOpEnd)
            return null

        // now we want to find out where the start and end is actually in relation to our sequence
        // reason is that the match is actually in relation to the aligned ops, but not to the
        // sequence we input


        // now we go backwards and work out the match sequence
        var numMatch = 0
        var seqIndex = -1
        var matchStart = -1
        var matchEnd = 0
        for (i in alignOps.indices)
        {
            val alignOp = alignOps[i]

            if (alignOp != SequenceAligner.AlignOp.DELETION)
            {
                // we need to keep track of where we are at the input sequence
                ++seqIndex
            }

            if (alignOp == SequenceAligner.AlignOp.MATCH && i >= alignOpStart && i < alignOpEnd)
            {
                ++numMatch

                if (matchStart == -1)
                {
                    matchStart = seqIndex
                }
                matchEnd = seqIndex + 1
            }
        }
        val matchSegmentLength = matchEnd - matchStart
        val matchRatio = numMatch.toDouble() / matchSegmentLength
        val matchSeq = seq.substring(matchStart, matchEnd)

        LOGGER.trace("seq({}) matchSegmentLength({}) numMatch({}) matchSegment({}) ratio({}) threshold({})",
            seq, matchSegmentLength, numMatch, matchSeq, matchRatio, matchThreshold)
        return TelomereMatch(matchStart, matchEnd, numMatch, matchSeq)
    }

    private fun generateTelomereTemplate(length: Int, gRich: Boolean): String
    {
        val canonicalHexamer = if (gRich) TealConstants.CANONICAL_TELOMERE_SEQ else TealConstants.CANONICAL_TELOMERE_SEQ_REV
        val b = StringBuilder(length)
        val numFullHexamers = length / canonicalHexamer.length
        for (i in 0 until numFullHexamers)
        {
            b.append(canonicalHexamer)
        }
        val residualLength = length % canonicalHexamer.length
        if (residualLength > 0)
        {
            b.append(canonicalHexamer, 0, residualLength)
        }
        return b.toString()
    }

    fun matchesGTelomere(seq: String, telomereMatchThreshold: Double, minTelomereMatchLength: Int): Boolean
    {
        return matchesTelomere(seq, telomereMatchThreshold, minTelomereMatchLength, gRich = true)
    }

    fun matchesCTelomere(seq: String, telomereMatchThreshold: Double, minTelomereMatchLength: Int): Boolean
    {
        return matchesTelomere(seq, telomereMatchThreshold, minTelomereMatchLength, gRich = false)
    }

    fun matchesTelomere(seq: String, telomereMatchThreshold: Double, minTelomereMatchLength: Int, gRich: Boolean): Boolean
    {
        if (!seq.contains(if (gRich) "GGG" else "CCC"))
        {
            return false
        }

        if (seq.length < minTelomereMatchLength)
        {
            // if is less than the length we require full match
            return findTelomereSegmentImpl(seq, 1.0, gRich) != null
        }

        val telomereMatch = findTelomereSegmentImpl(seq, telomereMatchThreshold, gRich)
        return telomereMatch != null && telomereMatch.numMatchedBases >= minTelomereMatchLength
    }
}
