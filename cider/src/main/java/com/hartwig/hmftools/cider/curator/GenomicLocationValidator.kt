package com.hartwig.hmftools.cider.curator

import com.hartwig.hmftools.cider.genes.GenomicLocation
import com.hartwig.hmftools.common.genome.region.Strand
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import htsjdk.samtools.util.SequenceUtil.reverseComplement
import org.apache.logging.log4j.LogManager

class GenomicLocationValidator(val refGenome: IndexedFastaSequenceFile)
{
    fun correctGenomicLocation(seq: String, genomicLocation: GenomicLocation, minBasesToCompare: Int, maxBasesToCompare: Int)
            : GenomicLocation
    {
        // get the ref sequence
        val refSeq = queryRefSequence(refGenome, genomicLocation)
        return correctGenomicLocation(seq, refSeq, genomicLocation, minBasesToCompare, maxBasesToCompare)
    }

    fun correctGenomicLocation(seq: String, refSeq: String, genomicLocation: GenomicLocation, minBasesToCompare: Int, maxBasesToCompare: Int)
            : GenomicLocation
    {
        // the sequence from the two files are not exact matches, so we need to correct for it
        // note that we only try to find the first section of the sequence in ref seq
        // this is to ensure it works even if there are some bases inserted / deleted somewhere
        val (startShift, numStartMismatch) = calcPositionMismatch(seq, refSeq, minBasesToCompare, maxBasesToCompare)
        var (endShift, numEndMismatch) = calcPositionMismatch(seq.reversed(), refSeq.reversed(), minBasesToCompare, maxBasesToCompare)
        endShift = -endShift

        if (numStartMismatch >= 10 || numEndMismatch >= 10)
        {
            // the sequences between ref and the location we got from IMGT is not a full match. This would be ok as long as
            // the anchor location is a good match
            sLogger.warn("mismatch count({}) non matching sequence: seq({}) refSeq({})", numStartMismatch, seq, refSeq)
        }

        if (startShift == 0 && seq.length == refSeq.length)
            return genomicLocation

        val correctedGeneLoc = if (genomicLocation.strand == Strand.FORWARD)
            genomicLocation.copy(posStart = genomicLocation.posStart + startShift,
                posEnd = genomicLocation.posEnd + endShift)
        else
            genomicLocation.copy(
                posStart = genomicLocation.posStart - endShift,
                posEnd = genomicLocation.posEnd - startShift)

        sLogger.info("seq({}) refSeq({}) refLocation({}) shift({}, {}) corrected({})",
            seq, refSeq, genomicLocation, startShift, endShift, correctedGeneLoc)

        return correctedGeneLoc
    }

    fun calcPositionMismatch(seq: String, refSeq: String, minBasesToCompare: Int, maxBasesToCompare: Int): Pair<Int, Int>
    {
        val minBasesToCompare = minOf(seq.length, refSeq.length, minBasesToCompare)

        var bestShift = 0
        var lowestNumMismatch = Int.MAX_VALUE

        val minShift = -seq.length + minBasesToCompare
        val maxShift = refSeq.length - minBasesToCompare

        // we shift along to find the best match
        for (i in minShift..maxShift)
        {
            var j = Math.max(0, i)
            var numMismatch = 0
            var numCompare = 0
            while ((j - i) < seq.length && j < refSeq.length && numCompare < maxBasesToCompare)
            {
                if (seq[j - i] != refSeq[j])
                {
                    ++numMismatch
                }
                ++numCompare
                ++j
            }

            if (numCompare >= minBasesToCompare && numMismatch < lowestNumMismatch)
            {
                lowestNumMismatch = numMismatch
                bestShift = i
            }
        }

        return Pair(bestShift, lowestNumMismatch)
    }

    // validate sequence against the ref genome file to make sure we got it right
    fun validateAgainstRefGenome(seq: String, genomicLocation: GenomicLocation): Boolean
    {
        val refGenomeSeq = queryRefSequence(refGenome, genomicLocation)

        if (refGenomeSeq.length != seq.length)
        {
            sLogger.warn("validation failed: seq({}) and ref genome seq({} of {}) length mismatch", seq, refGenomeSeq, genomicLocation)
            return false
        }

        // do not allow more than 2 bases difference
        var numDiff: Int = 0
        for (i in seq.indices)
        {
            if (seq[i] != refGenomeSeq[i])
                ++numDiff
        }

        if (numDiff > 6)
        {
            sLogger.error(
                "validation failed: seq({}) and ref genome seq({} of {}) sequence mismatch({}) > 6",
                seq,
                refGenomeSeq,
                genomicLocation,
                numDiff
            )
            return false
        }
        return true
    }

    fun queryRefSequence(refGenome: IndexedFastaSequenceFile, genomicLocation: GenomicLocation): String
    {
        var chromosome = genomicLocation.chromosome

        if (!refGenome.index.hasIndexEntry(chromosome))
        {
            // maybe need to try removing chr
            chromosome = chromosome.replace("chr", "")
        }

        val forwardSeq = refGenome.getSubsequenceAt(
            chromosome,
            genomicLocation.posStart.toLong(), genomicLocation.posEnd.toLong()
        ).baseString
        return if (genomicLocation.strand == Strand.FORWARD)
            forwardSeq
        else
            reverseComplement(forwardSeq)
    }

    companion object
    {
        private val sLogger = LogManager.getLogger(GenomicLocationValidator::class.java)
    }
}