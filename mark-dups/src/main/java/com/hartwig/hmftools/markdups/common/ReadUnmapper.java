package com.hartwig.hmftools.markdups.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.loadChrBaseRegions;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CHROMOSOME_INDEX;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.UNMAP_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.mateUnmapped;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;

import htsjdk.samtools.SAMRecord;

public class ReadUnmapper
{
    private final Map<String,List<BaseRegion>> mChrLocationsMap; // keyed by chromosome start
    private boolean mEnabled;
    private final UnmapStats mStats;

    private static final double MIN_OVERLAP_PERC = 0.9;

    public ReadUnmapper(final String filename)
    {
        this(loadChrBaseRegions(filename));

        if(mEnabled)
        {
            MD_LOGGER.info("loaded {} unmapped regions from {}",
                    mChrLocationsMap.values().stream().mapToInt(x -> x.size()).sum(), filename);
        }
    }

    public ReadUnmapper(final Map<String,List<BaseRegion>> chrLocationsMap)
    {
        mChrLocationsMap = chrLocationsMap;
        mEnabled = mChrLocationsMap != null && !mChrLocationsMap.isEmpty();
        mStats = new UnmapStats();
    }

    public List<BaseRegion> getRegions(final String chromosome) { return mChrLocationsMap.get(chromosome); }
    public boolean enabled() { return mEnabled; }
    public UnmapStats stats() { return mStats; }

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
            unmapMate = mateReadInRegion(read);
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

    @VisibleForTesting
    public static boolean readInRegion(final SAMRecord read, final List<BaseRegion> regions)
    {
        int requiredBaseCount = (int)((read.getAlignmentEnd() - read.getAlignmentStart()) * MIN_OVERLAP_PERC);
        return containsReadPositions(read.getAlignmentStart(), read.getAlignmentEnd(), requiredBaseCount, regions);
    }

    @VisibleForTesting
    public boolean mateReadInRegion(final SAMRecord read)
    {
        List<BaseRegion> mateRegions = mChrLocationsMap.get(read.getMateReferenceName());

        if(mateRegions == null)
            return false;

        // coordinate check must be identical to how the mate checks itself, which requires knowledge of aligned bases
        if(read.hasAttribute(MATE_CIGAR_ATTRIBUTE))
        {
            int approxMateEnd = read.getMateAlignmentStart() + read.getReadBases().length - 1;

            return containsReadPositions(
                    read.getMateAlignmentStart(), approxMateEnd, read.getStringAttribute(MATE_CIGAR_ATTRIBUTE), mateRegions);
        }
        else
        {
            // TODO: this won't work for split reads
            int readBaseLength = read.getReadBases().length;

            return containsReadPositions(
                    read.getMateAlignmentStart(), read.getMateAlignmentStart() + readBaseLength - 1, readBaseLength, mateRegions);
        }
    }

    @VisibleForTesting
    public boolean supplementaryInRegion(final SAMRecord read)
    {
        SupplementaryReadData suppReadData = SupplementaryReadData.firstAlignmentFrom(read);

        if(suppReadData == null)
            return false;

        List<BaseRegion> suppRegions = mChrLocationsMap.get(suppReadData.Chromosome);

        if(suppRegions == null)
            return false;

        int approxSuppEnd = suppReadData.Position + read.getReadBases().length - 1;

        return containsReadPositions(suppReadData.Position, approxSuppEnd, suppReadData.Cigar, suppRegions);
    }

    private static boolean containsReadPositions(
            final int readStart, final int readEnd, final int requiredReadBaseCount, final List<BaseRegion> regions)
    {
        for(BaseRegion region : regions)
        {
            if(!positionsOverlap(readStart, readEnd, region.start(), region.end()))
                continue;

            if(positionsWithin(readStart, readEnd, region.start(), region.end()))
                return true;

            int overlapBases = min(region.end(), readEnd) - max(region.start(), readStart) + 1;

            int requiredRegionBaseCount = (int)(region.baseLength() * MIN_OVERLAP_PERC);

            if(overlapBases >= requiredReadBaseCount || overlapBases >= requiredRegionBaseCount)
                return true;
        }

        return false;
    }

    private boolean containsReadPositions(final int readStart, int approxReadEnd, final String cigar, final List<BaseRegion> regions)
    {
        // mate coords and supplementaries don't have alignment end so this needs to be determined from the cigar, but only
        // do this is if they are a possible overlap with a region

        int readEnd = -1;
        int requiredReadBaseCount = 0;

        for(BaseRegion region : regions)
        {
            if(!positionsOverlap(readStart, approxReadEnd, region.start(), region.end()))
                continue;

            // now test with precise read end alignment
            if(readEnd < 0)
            {
                readEnd = getReadEndFromCigar(readStart, cigar);
                requiredReadBaseCount = (int)((readEnd - readStart) * MIN_OVERLAP_PERC);
            }

            int overlapBases = min(region.end(), readEnd) - max(region.start(), readStart) + 1;

            int requiredRegionBaseCount = (int)(region.baseLength() * MIN_OVERLAP_PERC);

            if(overlapBases >= requiredReadBaseCount || overlapBases >= requiredRegionBaseCount)
                return true;
        }

        return false;
    }

    public int getReadEndFromCigar(final int readStart, final String cigarStr)
    {
        int currentPosition = readStart;
        int elementLength = 0;

        for(int i = 0; i < cigarStr.length(); ++i)
        {
            char c = cigarStr.charAt(i);
            boolean isAddItem = (c == 'D' || c == 'M' || c == 'N');

            if(isAddItem)
            {
                currentPosition += elementLength;
                elementLength = 0;
                continue;
            }

            int digit = c - '0';
            if (digit >= 0 && digit <= 9)
            {
                elementLength = elementLength * 10 + digit;
            }
            else
            {
                elementLength = 0;
            }
        }

        // always pointing to the start of the next element, so need to move back a base
        return currentPosition - 1;
    }

    private static void setUnmappedAttributes(final SAMRecord read)
    {
        read.setReadUnmappedFlag(true);
        read.setMappingQuality(0);

        // store the original mapping
        setUnmapCoordsAttribute(read, read.getReferenceName(), read.getAlignmentStart());
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

    private static final String UNMAPP_COORDS_DELIM = ":";

    public static String[] parseUnmappedCoords(final SAMRecord read)
    {
        return read.getStringAttribute(UNMAP_ATTRIBUTE).split(UNMAPP_COORDS_DELIM, 2);
    }

    private static void setUnmapCoordsAttribute(final SAMRecord read, final String chromosome, final int position)
    {
        read.setAttribute(UNMAP_ATTRIBUTE, chromosome + UNMAPP_COORDS_DELIM + position);
    }

    public static void unmapMateAlignment(final SAMRecord read, boolean mateUnmapped)
    {
        // set flag mate unmapped
        read.setMateUnmappedFlag(true);
        read.setProperPairFlag(false);

        // store the mate's original mapping
        setUnmapCoordsAttribute(read, read.getMateReferenceName(), read.getMateAlignmentStart());

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

    @VisibleForTesting
    public void addRegion(final String chromosome, final BaseRegion region)
    {
        List<BaseRegion> regions = mChrLocationsMap.get(chromosome);

        if(regions == null)
        {
            regions = Lists.newArrayList();
            mChrLocationsMap.put(chromosome, regions);
        }

        regions.add(region);

        mEnabled = true;
    }
}
