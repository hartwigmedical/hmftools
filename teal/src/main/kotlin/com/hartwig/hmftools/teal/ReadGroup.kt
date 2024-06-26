package com.hartwig.hmftools.teal

import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SAMTag
import com.hartwig.hmftools.common.region.ChrBaseRegion
import com.hartwig.hmftools.common.bam.SamRecordUtils
import org.apache.logging.log4j.Level
import org.apache.logging.log4j.LogManager
import java.util.ArrayList

class ReadGroup(val name: String)
{
    data class SupplementaryAlignment(
        val firstOfPair: Boolean,
        val chromosome: String,
        val negativeStrand: Boolean,
        val position: Int,
        val cigar: String)
    {
        fun isMatch(read: SAMRecord): Boolean
        {
            return SamRecordUtils.firstInPair(read) == firstOfPair &&
                read.alignmentStart == position &&
                read.readNegativeStrandFlag == negativeStrand &&
                read.referenceName == chromosome &&
                cigarEqual(read.cigarString, cigar)
        }

        companion object
        {
            // for the purpose of matching supplementary alignments, hard clip and soft clip are equivalent
            private fun cigarEqual(cigar1: String, cigar2: String?): Boolean
            {
                return cigar1.replace('H', 'S') == cigar2!!.replace('H', 'S')
            }
        }
    }

    val mutableReads: MutableList<SAMRecord> = ArrayList()
    val mutableSupplementaryReads: MutableList<SAMRecord> = ArrayList()

    val Reads: List<SAMRecord> get() { return mutableReads }
    val SupplementaryReads: List<SAMRecord> get() { return mutableSupplementaryReads }

    // we have to find it
    val firstOfPair: SAMRecord?
        get()
        {
            // we have to find it
            for (r in Reads)
            {
                if (SamRecordUtils.firstInPair(r))
                {
                    return r
                }
            }
            return null
        }

    // we have to find it
    val secondOfPair: SAMRecord?
        get()
        {
            // we have to find it
            for (r in Reads)
            {
                if (!SamRecordUtils.firstInPair(r))
                {
                    return r
                }
            }
            return null
        }

    // get all the reads including non supplementary
    // and supplementary reads
    val allReads: List<SAMRecord>
        get()
        {
            return Reads + SupplementaryReads
        }

    fun isComplete(logLevel: Level? = null): Boolean
    {
        if (Reads.isEmpty())
        {
            logger.log(logLevel, "Read is empty")
            return false
        }

        // we check several things
        if (Reads[0].readPairedFlag && Reads.size != 2)
        {
            // we havent got all the reads yet
            logger.log(logLevel, "{} missing mate pair", Reads[0])
            return false
        }

        // check for any of the supplementary reads
        for (read in Reads)
        {
            val saAttribute = read.getStringAttribute(SAMTag.SA.name)
            if (saAttribute != null)
            {
                for (sa in suppAlignmentPositions(SamRecordUtils.firstInPair(read), saAttribute)!!)
                {
                    // check if this supplementary read exists
                    if (SupplementaryReads.stream().noneMatch { r: SAMRecord -> sa.isMatch(r) })
                    {
                        logger.log(
                            logLevel, "{} Missing supplementary read: aligned to {}:{}", Reads[0],
                            sa.chromosome, sa.position)
                        return false
                    }
                }
            }
        }
        return true
    }

    operator fun contains(read: SAMRecord): Boolean
    {
        if (read.readName != name)
        {
            return false
        }
        val listToLook: List<SAMRecord>
        listToLook = if (read.isSecondaryOrSupplementary)
        {
            SupplementaryReads
        }
        else
        {
            Reads
        }

        // we only check chromosome and alignment start
        return listToLook.stream()
            .anyMatch { x: SAMRecord ->
                    SamRecordUtils.firstInPair(x) == SamRecordUtils.firstInPair(read) &&
                    x.alignmentStart == read.alignmentStart &&
                    x.readNegativeStrandFlag == read.readNegativeStrandFlag &&
                    x.referenceName == read.referenceName &&
                    x.cigar == read.cigar }
    }

    //public boolean isDuplicate() { return Reads.stream().anyMatch(x -> x.isDuplicate()); }
    override fun toString(): String
    {
        return String.format("%s reads(%d) complete(%s)", name, Reads.size, isComplete())
    }

    fun invariant(): Boolean
    {
        // check to make sure several things:
        // 1. same record cannot appear more than once
        // 2. Reads cannot contain supplementary
        // 3. SupplementaryReads must only contain supplementary
        if (Reads.size > 2)
        {
            return false
        }
        if (!Reads.stream().allMatch { x: SAMRecord -> x.readName == name })
        {
            return false
        }
        if (!SupplementaryReads.stream().allMatch { x: SAMRecord -> x.readName == name })
        {
            return false
        }
        return if (Reads.stream().anyMatch { obj: SAMRecord -> obj.supplementaryAlignmentFlag })
        {
            false
        }
        else SupplementaryReads.stream().allMatch { x: SAMRecord -> x.supplementaryAlignmentFlag }
    }

    fun findMissingReadBaseRegions(): List<ChrBaseRegion>
    {
        val baseRegions: MutableList<ChrBaseRegion> = ArrayList()
        assert(invariant())
        if (Reads.size == 1)
        {
            val read = Reads[0]
            if (read.readPairedFlag && read.mateReferenceIndex != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX)
            {
                // note that even if the mate unmapped flag is set, we still might not be there
                // unmapped read of a mapped read pair would show up with the mate position
                if (read.mateReferenceName.isEmpty() || read.mateAlignmentStart <= 0)
                {
                    // this shouldn't happen
                    logger.warn(
                        "read({}) invalid mate reference, mate ref index({})",
                        read, read.mateReferenceIndex)
                }
                else
                {
                    baseRegions.add(ChrBaseRegion(read.mateReferenceName, read.mateAlignmentStart, read.mateAlignmentStart))
                    logger.trace(
                        "{} missing read mate: aligned to {}:{}", Reads[0],
                        read.mateReferenceName, read.mateAlignmentStart)
                }
            }
        }
        for (read in Reads)
        {
            val saAttribute = read.getStringAttribute(SAMTag.SA.name)
            if (saAttribute != null)
            {
                for (sa in suppAlignmentPositions(SamRecordUtils.firstInPair(read), saAttribute)!!)
                {
                    // check if this supplementary read exists
                    if (SupplementaryReads.stream().noneMatch { r: SAMRecord -> sa.isMatch(r) })
                    {
                        baseRegions.add(ChrBaseRegion(sa.chromosome, sa.position, sa.position))
                        logger.trace(
                            "{} Missing supplementary read: aligned to {}:{}", read,
                            sa.chromosome, sa.position)
                    }
                }
            }
        }
        return baseRegions
    }

    fun acceptRead(record: SAMRecord): Boolean
    {
        if (contains(record)) return false
        if (record.isSecondaryOrSupplementary)
        {
            // bam files generated by some GATK IndelRealigner
            // uses secondary instead of supplementary alignment
            if (record.isSecondaryAlignment)
            {
                record.supplementaryAlignmentFlag = true
                record.isSecondaryAlignment = false
            }
            mutableSupplementaryReads.add(record)
        }
        else
        {
            mutableReads.add(record)
        }
        return true
    }

    companion object
    {
        private val logger = LogManager.getLogger(ReadGroup::class.java)
        
        const val SUPP_ALIGNMENT_DELIM = ","
        @JvmStatic
        fun suppAlignmentPositions(firstOfPair: Boolean, suppAlignment: String?): List<SupplementaryAlignment>?
        {
            if (suppAlignment == null) return null
            val suppAignPos: MutableList<SupplementaryAlignment> = ArrayList()
            val supplementaryItems = suppAlignment.split(";").toTypedArray()
            for (si in supplementaryItems)
            {
                val items = si.split(SUPP_ALIGNMENT_DELIM).toTypedArray()
                if (items.size < 5) continue

                // supplementary(SA) string attribute looks like
                // SA:Z:(rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+
                // 20,61647163,+,99M52S,0,1;11,70524575,+,95S30M26S,0,0;
                // the first word is the chromosome, the second is the alignment start
                val sa = SupplementaryAlignment(
                    firstOfPair = firstOfPair,
                    chromosome = items[0],
                    position = items[1].toInt(),
                    negativeStrand = items[2] == "-",
                    cigar = items[3]
                )
                suppAignPos.add(sa)
            }
            return suppAignPos
        }
    }
}