package com.hartwig.hmftools.bamtools.btofmc;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public interface RegionTask
{
    RegionTask UNMAPPED_READS = new UnmappedReads();

    static RegionTask createChrRegion(ChrBaseRegion region)
    {
        return new ChrRegion(region);
    }

    void slice(final SamReader samReader, final BamSlicer bamSlicer, final Consumer<SAMRecord> consumer);

    boolean isAlignmentStartWithinRegion(final SAMRecord read);

    class UnmappedReads implements RegionTask
    {
        private UnmappedReads()
        {
        }

        public void slice(final SamReader samReader, final BamSlicer bamSlicer, final Consumer<SAMRecord> consumer)
        {
            bamSlicer.queryUnmapped(samReader, consumer);
        }

        public boolean isAlignmentStartWithinRegion(final SAMRecord read)
        {
            return true;
        }

        @Override
        public String toString()
        {
            return "UnmappedReads";
        }
    }

    class ChrRegion implements RegionTask
    {
        private final ChrBaseRegion mChrRegion;

        public ChrRegion(final ChrBaseRegion chrRegion)
        {
            mChrRegion = chrRegion;
        }

        public void slice(final SamReader samReader, final BamSlicer bamSlicer, final Consumer<SAMRecord> consumer)
        {
            bamSlicer.slice(samReader, mChrRegion, consumer);
        }

        public boolean isAlignmentStartWithinRegion(final SAMRecord read)
        {
            return mChrRegion.containsPosition(read.getAlignmentStart());
        }

        @Override
        public String toString()
        {
            return mChrRegion.toString();
        }
    }
}
