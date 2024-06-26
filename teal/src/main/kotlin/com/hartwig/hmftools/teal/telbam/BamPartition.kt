package com.hartwig.hmftools.teal.telbam

import com.hartwig.hmftools.common.region.ChrBaseRegion
import htsjdk.samtools.SAMFileHeader
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SAMRecordIterator
import htsjdk.samtools.SamReader
import org.apache.commons.lang3.Validate

// abstraction of a bam partition
interface BamPartition
{
    override fun toString(): String
    fun iterator(samReader: SamReader): SAMRecordIterator

    // partition by genomic region. It only return reads where the alignment start is
    // inside the region
    class BamPartitionByRegion(private val mBaseRegion: ChrBaseRegion) : BamPartition
    {
        override fun toString(): String
        {
            return String.format("%s:%,d-%,d", mBaseRegion.chromosome(), mBaseRegion.start(), mBaseRegion.end())
        }

        override fun iterator(samReader: SamReader): SAMRecordIterator
        {
            val overlappingItr = samReader.query(mBaseRegion.chromosome(), mBaseRegion.start(), mBaseRegion.end(), false)

            // wrap the iterator to skip over records where the start alignment is not inside the region
            // this is to avoid processing a bam record twice if it overlaps two regions
            return object : SAMRecordIterator
            {
                var nextRead: SAMRecord? = null
                override fun assertSorted(sortOrder: SAMFileHeader.SortOrder): SAMRecordIterator
                {
                    return overlappingItr.assertSorted(sortOrder)
                }

                override fun close()
                {
                    overlappingItr.close()
                }

                override fun remove()
                {
                    throw UnsupportedOperationException("remove")
                }

                override fun hasNext(): Boolean
                {
                    if (nextRead != null) return true
                    while (overlappingItr.hasNext())
                    {
                        val r = overlappingItr.next()
                        if (mBaseRegion.containsPosition(r.alignmentStart))
                        {
                            nextRead = r
                            break
                        }
                    }
                    return nextRead != null
                }

                override fun next(): SAMRecord
                {
                    return if (nextRead != null || hasNext())
                    {
                        val r = nextRead
                        nextRead = null
                        assert(r != null)
                        Validate.isTrue(mBaseRegion.containsPosition(r!!.alignmentStart))
                        r
                    } else throw NoSuchElementException()
                }
            }
        }
    }

    // all unmapped reads
    class BamPartitionUnmapped : BamPartition
    {
        override fun toString(): String
        {
            return "unmapped"
        }

        override fun iterator(samReader: SamReader): SAMRecordIterator
        {
            return samReader.queryUnmapped()
        }
    }

    // whole bam as one partition
    class BamPartitionWholeBam : BamPartition
    {
        override fun toString(): String
        {
            return "wholeBam"
        }

        override fun iterator(samReader: SamReader): SAMRecordIterator
        {
            return samReader.iterator()
        }
    }

    companion object
    {
        // a bam partition of reads that start in region
        fun ofRegion(baseRegion: ChrBaseRegion): BamPartition
        {
            return BamPartitionByRegion(baseRegion)
        }

        fun ofUnmapped(): BamPartition
        {
            return BamPartitionUnmapped()
        }

        fun ofWholeBam(): BamPartition
        {
            return BamPartitionWholeBam()
        }
    }
}