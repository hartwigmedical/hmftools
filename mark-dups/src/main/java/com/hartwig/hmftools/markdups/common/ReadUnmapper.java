package com.hartwig.hmftools.markdups.common;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.loadChrBaseRegions;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CHROMOSOME_INDEX;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.UNMAP_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.mateUnmapped;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;

import htsjdk.samtools.SAMRecord;

public class ReadUnmapper
{
    private final Map<String,List<BaseRegion>> mChrLocationsMap; // keyed by chromosome start
    private final int mReadLength;
    private boolean mEnabled;
    private final UnmapStats mStats;

    public ReadUnmapper(final String filename, final int readLength)
    {
        this(loadChrBaseRegions(filename), readLength);

        if(mEnabled)
        {
            MD_LOGGER.info("loaded {} unmapped regions from {}",
                    mChrLocationsMap.values().stream().mapToInt(x -> x.size()).sum(), filename);
        }
    }

    public ReadUnmapper(final Map<String,List<BaseRegion>> chrLocationsMap, final int readLength)
    {
        mReadLength = readLength;
        mChrLocationsMap = chrLocationsMap;
        mEnabled = mChrLocationsMap != null && !mChrLocationsMap.isEmpty();
        mStats = new UnmapStats();
    }

    public List<BaseRegion> getRegions(final String chromosome) { return mChrLocationsMap.get(chromosome); }
    public boolean enabled() { return mEnabled; }
    public UnmapStats stats() { return mStats; }

    public static boolean readInRegion(final SAMRecord read, final List<BaseRegion> regions)
    {
        return containsReadPositions(read.getAlignmentStart(), read.getAlignmentEnd(), regions);
    }

    private static boolean containsReadPositions(final int readStart, final int readEnd, final List<BaseRegion> regions)
    {
        return regions.stream().anyMatch(x -> positionsWithin(readStart, readEnd, x.start(), x.end()));
    }

    public boolean mateReadInRegion(final SAMRecord read)
    {
        List<BaseRegion> mateRegions = mChrLocationsMap.get(read.getMateReferenceName());

        if(mateRegions == null)
            return false;

        return containsReadPositions(read.getMateAlignmentStart(), read.getMateAlignmentStart() + mReadLength - 1, mateRegions);
    }

    public boolean supplementaryInRegion(final SAMRecord read)
    {
        SupplementaryReadData suppReadData = SupplementaryReadData.from(read);

        if(suppReadData == null)
            return false;

        List<BaseRegion> suppRegions = mChrLocationsMap.get(suppReadData.Chromosome);

        if(suppRegions == null)
            return false;

        return containsReadPositions(suppReadData.Position, suppReadData.Position + mReadLength - 1, suppRegions);
    }

    public boolean checkTransformRead(final SAMRecord read, final List<BaseRegion> readRegions)
    {
        if(!mEnabled || readRegions == null)
            return false;

        // first check the read's alignment itself
        boolean readUnmapped = read.getReadUnmappedFlag();
        boolean mateUnmapped = mateUnmapped(read);

        if(readUnmapped && mateUnmapped)
            return false; // nothing to do and there won't be a supplementary

        boolean unmapRead = !readUnmapped && readInRegion(read, readRegions);

        if(unmapRead && read.getSupplementaryAlignmentFlag())
        {
            // these will be dropped from the BAM
            setUnmappedAttributes(read);
            mStats.SupplementaryCount.incrementAndGet();
            return true;
        }

        boolean unmapMate = false;

        if(!mateUnmapped)
        {
            List<BaseRegion> mateRegions = mChrLocationsMap.get(read.getMateReferenceName());
            unmapMate = mateRegions != null && mateReadInRegion(read);
        }

        if(unmapRead)
        {
            unmapReadAlignment(read, mateUnmapped || unmapMate);
        }

        if(unmapMate)
        {
            unmapMateAlignment(read, readUnmapped || unmapRead);
        }

        if((readUnmapped || unmapRead) && (mateUnmapped || unmapMate))
        {
            mStats.UnmappedCount.incrementAndGet();
        }
        else if(unmapRead)
        {
            mStats.ReadCount.incrementAndGet();
        }
        else if(unmapMate)
        {
            mStats.MateCount.incrementAndGet();
        }

        // check supplementary
        boolean unmapSuppAlignment = false;
        if(!unmapRead && supplementaryInRegion(read))
        {
            unmapSupplementary(read);
            mStats.SuppAlignmentCount.incrementAndGet();
            unmapSuppAlignment = true;
        }

        return unmapRead || unmapMate || unmapSuppAlignment;
    }

    private static void setUnmappedAttributes(final SAMRecord read)
    {
        read.setReadUnmappedFlag(true);
        read.setMappingQuality(0);
        read.setAttribute(UNMAP_ATTRIBUTE, 1);
    }

    public static void unmapReadAlignment(final SAMRecord read, boolean mateUnmapped)
    {
        setUnmappedAttributes(read);

        read.setProperPairFlag(false);

        // clear insert size
        read.setInferredInsertSize(0);

        // clear reference index, reference name and alignment start
        if(mateUnmapped)
        {
            // set both to unknown
            read.setAlignmentStart(0);
            read.setReferenceIndex(NO_CHROMOSOME_INDEX);
            read.setReferenceName(NO_CHROMOSOME_NAME);
        }
        else
        {
            // set to primary's details
            read.setAlignmentStart(read.getMateAlignmentStart());
            read.setReferenceIndex(read.getMateReferenceIndex());
            read.setReferenceName(read.getMateReferenceName());
        }

        // clear cigar if present
        read.setCigarString(NO_CIGAR);

        // clear the supplementary data too
        unmapSupplementary(read);
    }

    public static void unmapMateAlignment(final SAMRecord read, boolean mateUnmapped)
    {
        // set flag mate unmapped
        read.setMateUnmappedFlag(true);
        read.setProperPairFlag(false);

        // clear insert size
        read.setInferredInsertSize(0);

        // set mate reference index, reference name and alignment start
        if(mateUnmapped)
        {
            // set both to unknown
            read.setMateAlignmentStart(0);
            read.setMateReferenceIndex(NO_CHROMOSOME_INDEX);
            read.setMateReferenceName(NO_CHROMOSOME_NAME);
        }
        else
        {
            // set to primary's details
            read.setMateAlignmentStart(read.getAlignmentStart());
            read.setMateReferenceIndex(read.getReferenceIndex());
            read.setMateReferenceName(read.getReferenceName());
        }

        // clear mate cigar if present
        if(read.hasAttribute(MATE_CIGAR_ATTRIBUTE))
        {
            read.setAttribute(MATE_CIGAR_ATTRIBUTE, null);
        }
    }

    public static void unmapSupplementary(final SAMRecord read)
    {
        // remove supplementary alignment details from attributes
        if(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
        {
            read.setAttribute(SUPPLEMENTARY_ATTRIBUTE, null);
        }
    }
}
