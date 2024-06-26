package com.hartwig.hmftools.bamtools.compare;

import java.util.NoSuchElementException;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.commons.lang3.Validate;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

// abstraction of a bam partition
public interface BamPartition
{
    String toString();
    SAMRecordIterator iterator(SamReader samReader);

    // a bam partition of reads that start in region
    static BamPartition ofRegion(final ChrBaseRegion baseRegion)
    {
        return new BamPartitionByRegion(baseRegion);
    }

    static BamPartition ofUnmapped()
    {
        return new BamPartitionUnmapped();
    }

    static BamPartition ofWholeBam()
    {
        return new BamPartitionWholeBam();
    }

    // partition by genomic region. It only return reads where the alignment start is
    // inside the region
    class BamPartitionByRegion implements BamPartition
    {
        private final ChrBaseRegion mBaseRegion;

        public BamPartitionByRegion(final ChrBaseRegion baseRegion)
        {
            mBaseRegion = baseRegion;
        }

        @Override
        public String toString()
        {
            return String.format("%s:%,d-%,d", mBaseRegion.chromosome(), mBaseRegion.start(), mBaseRegion.end());
        }

        @Override
        public SAMRecordIterator iterator(SamReader samReader)
        {
            SAMRecordIterator overlappingItr = samReader.query(mBaseRegion.chromosome(),
                    mBaseRegion.start(), mBaseRegion.end(), false);

            // wrap the iterator to skip over records where the start alignment is not inside the region
            // this is to avoid processing a bam record twice if it overlaps two regions
            return new SAMRecordIterator()
            {
                SAMRecord nextRead = null;
                @Override
                public SAMRecordIterator assertSorted(final SAMFileHeader.SortOrder sortOrder)
                {
                    return overlappingItr.assertSorted(sortOrder);
                }
                @Override
                public void close() { overlappingItr.close(); }
                @Override
                public boolean hasNext()
                {
                    if(nextRead != null)
                        return true;
                    while(overlappingItr.hasNext())
                    {
                        SAMRecord r = overlappingItr.next();
                        if(mBaseRegion.containsPosition(r.getAlignmentStart()))
                        {
                            nextRead = r;
                            break;
                        }
                    }
                    return nextRead != null;
                }
                @Override
                public SAMRecord next()
                {
                    if(nextRead != null || hasNext())
                    {
                        SAMRecord r = nextRead;
                        nextRead = null;
                        assert r != null;
                        Validate.isTrue(mBaseRegion.containsPosition(r.getAlignmentStart()));
                        return r;
                    }
                    else
                        throw new NoSuchElementException();
                }
            };
        }
    }

    // all unmapped reads
    class BamPartitionUnmapped implements BamPartition
    {
        @Override
        public String toString()
        {
            return "unmapped";
        }

        @Override
        public SAMRecordIterator iterator(SamReader samReader)
        {
            return samReader.queryUnmapped();
        }
    }

    // whole bam as one partition
    class BamPartitionWholeBam implements BamPartition
    {
        @Override
        public String toString()
        {
            return "wholeBam";
        }

        @Override
        public SAMRecordIterator iterator(SamReader samReader)
        {
            return samReader.iterator();
        }
    }
}
