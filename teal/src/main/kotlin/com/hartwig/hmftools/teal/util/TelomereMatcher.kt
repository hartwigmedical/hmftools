package com.hartwig.hmftools.teal.util

import com.hartwig.hmftools.common.aligner.*

import com.hartwig.hmftools.teal.TealConstants
import org.apache.logging.log4j.LogManager

// helper class to match a sequence to telomere
// and give a similarity rate
object TelomereMatcher
{
    private val LOGGER = LogManager.getLogger(javaClass)

    // we use same score for mismatch and gap, tailored for aligning to TTAGGG template
    private val sAligner = LocalSequenceAligner(1, -1, -1, -1)

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
        val alignment: LocalSequenceAligner.Alignment = sAligner.alignSequence(seq, generateTelomereTemplate(telomereTemplateLength, gRich))

        val numMatch = alignment.operators.count({ op -> op == AlignmentOperator.MATCH })

        // at least need to match 12 (2 x TTAGGG)
        if (numMatch < TealConstants.MIN_TELOMERE_MATCH_BASES)
            return null

        val matchRatio = numMatch / alignment.firstSequenceAlignLength
        val matchSeq = seq.substring(alignment.firstSequenceAlignStart, alignment.firstSequenceAlignEnd)

        LOGGER.trace("seq({}) matchSegmentLength({}) numMatch({}) matchSegment({}) ratio({}) threshold({})",
            seq, alignment.firstSequenceAlignLength, numMatch, matchSeq, matchRatio, matchThreshold)
        return TelomereMatch(alignment.firstSequenceAlignStart, alignment.firstSequenceAlignEnd, numMatch, matchSeq)
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
