package com.hartwig.hmftools.cider.curator

import com.hartwig.hmftools.cider.genes.GenomicLocation
import com.hartwig.hmftools.common.genome.region.Strand
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import htsjdk.samtools.util.SequenceUtil.reverseComplement
import org.apache.logging.log4j.LogManager

class GenomicLocationValidator(val refGenome: IndexedFastaSequenceFile)
{
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
        var contig = genomicLocation.contig.name

        if (!refGenome.index.hasIndexEntry(contig))
        {
            // maybe need to try removing chr
            contig = contig.replace("chr", "")
        }

        val forwardSeq = refGenome.getSubsequenceAt(
            contig,
            genomicLocation.position.start().toLong(), genomicLocation.position.end().toLong()
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