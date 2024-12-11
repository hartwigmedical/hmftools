package com.hartwig.hmftools.redux.unmap;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.BaseRegion.binarySearch;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_INDEX;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UNMAP_ATTRIBUTE;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.common.Constants.UNMAP_CHIMERIC_FRAGMENT_LENGTH_MAX;
import static com.hartwig.hmftools.redux.common.Constants.UNMAP_MAX_NON_OVERLAPPING_BASES;
import static com.hartwig.hmftools.redux.common.Constants.UNMAP_MIN_HIGH_DEPTH;
import static com.hartwig.hmftools.redux.common.Constants.UNMAP_MIN_SOFT_CLIP;

import java.util.List;
import java.util.Map;

import javax.annotation.Nullable;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.region.HighDepthRegion;
import com.hartwig.hmftools.common.region.UnmappedRegions;

import htsjdk.samtools.SAMRecord;

public class ReadUnmapper
{
    private final Map<String,List<HighDepthRegion>> mChrLocationsMap; // keyed by chromosome start
    private boolean mEnabled;
    private final UnmapStats mStats;

    public ReadUnmapper(final String filename)
    {
        this(loadUnmapRegions(filename));

        if(mEnabled)
        {
            RD_LOGGER.info("loaded {} unmapping regions from {}",
                    mChrLocationsMap.values().stream().mapToInt(x -> x.size()).sum(), filename);
        }
    }

    public ReadUnmapper(final Map<String,List<HighDepthRegion>> chrLocationsMap)
    {
        mChrLocationsMap = chrLocationsMap;
        mEnabled = mChrLocationsMap != null && !mChrLocationsMap.isEmpty();
        mStats = new UnmapStats();
    }

    public List<HighDepthRegion> getRegions(final String chromosome)
    {
        return mChrLocationsMap.get(chromosome);
    }
    public Map<String,List<HighDepthRegion>> getAllRegions() { return mChrLocationsMap; }

    public boolean enabled() { return mEnabled; }

    public UnmapStats stats() { return mStats; }

    public boolean checkTransformRead(final SAMRecord read, final UnmapRegionState regionState)
    {
        if(!mEnabled || regionState == null)
            return false;

        /* Criteria for unmapping a read:
            - falls within a region of high depth
            - discordant - INV, BND, one read unmapped or fragment length > 1000 (only for paired reads)
            - soft-clip bases > 20
            - non-human chromosome for now

           Scenarios & logic:
           - both read and mate are already unmapped, then nothing to do
           - check the read's coords, its mate's coords and any supplementary alignment coords vs the loaded unmapping regions
           - for any overlap with a non-high-depth region, additionally check discordant and soft-clip bases as above
           -
           - supplementaries - unmap if their primary or another associated supplementary will be, or if they need to be
         */

        if(read.isSecondaryAlignment())
            return checkSecondaryRead(read, regionState);
        else if(read.getReadPairedFlag())
            return checkTransformPairedRead(read, regionState);
        else
            return checkTransformUnpairedRead(read, regionState);
    }

    private boolean checkTransformPairedRead(final SAMRecord read, final UnmapRegionState regionState)
    {
        // first check the read's alignment itself
        boolean readUnmapped = read.getReadUnmappedFlag();
        boolean mateUnmapped = read.getMateUnmappedFlag();

        if(readUnmapped && mateUnmapped)
            return false; // nothing to do and there won't be a supplementary

        boolean isSupplementary = read.getSupplementaryAlignmentFlag();

        boolean unmapRead = false;
        boolean unmapMate = false;

        RegionMatchType readRegionType = !readUnmapped ? findMaxDepthRegionOverlap(
                read.getAlignmentStart(), read.getAlignmentEnd(), regionState.PartitionRegions, regionState, true) : null;

        UnmapReason unmapReason = null;

        if(!readUnmapped)
        {
            unmapReason = checkUnmapRead(read, readRegionType);
            unmapRead = unmapReason != UnmapReason.NONE;
        }

        RegionMatchType mateRegionType = !mateUnmapped ? mateMaxDepthRegionOverlap(read, regionState) : null;

        if(!mateUnmapped)
            unmapMate =  checkUnmapMate(read, mateRegionType);

        // handle scenario where read wouldn't be unmapped unless mate becomes unmapped
        if(!readUnmapped && unmapReason == UnmapReason.NONE && readRegionType != RegionMatchType.NONE && unmapMate)
        {
            unmapReason = UnmapReason.CHIMERIC;
            unmapRead = true;
        }

        if(!mateUnmapped && !unmapMate && mateRegionType != RegionMatchType.NONE && unmapRead)
        {
            unmapMate = true;
        }

        if(!readUnmapped)
        {
            if(unmapReason != UnmapReason.NONE)
            {
                switch(unmapReason)
                {
                    case HIGH_DEPTH:
                        mStats.HighDepthCount.incrementAndGet();
                        break;

                    case SOFT_CLIP:
                        mStats.LongSoftClipCount.incrementAndGet();
                        break;

                    case CHIMERIC:
                        mStats.ChimericCount.incrementAndGet();
                        break;

                    default:
                        break;
                }
            }
            else if(isSupplementary && !unmapRead)
            {
                // supplementaries are unmapped if their primary, or an associated supplementary, will be unmapped
                unmapRead = checkUnmapSupplementaryRead(read);
            }

            if(unmapRead && isSupplementary)
            {
                // these will be dropped from the BAM
                setUnmappedAttributes(read);
                mStats.SupplementaryCount.incrementAndGet();
                return true;
            }
        }

        if(unmapRead)
            unmapReadAlignment(read, mateUnmapped, unmapMate);

        if(unmapMate)
        {
            // Note that for supplementary reads, at this point readUnmapped will be false, since a supplementary must be mapped. Also,
            // unmapRead is false, else we would have returned already. Therefore, this will only modify the mate properties and attributes.
            unmapMateAlignment(read, readUnmapped, unmapRead);
        }

        if((readUnmapped || unmapRead) && (mateUnmapped || unmapMate))
            mStats.UnmappedCount.incrementAndGet();
        else if(unmapRead)
            mStats.ReadCount.incrementAndGet();
        else if(unmapMate)
            mStats.MateCount.incrementAndGet();

        // nothing more to do for a supplementary whose primary or mate has been unmapped, or if it was unmapped
        if(isSupplementary)
            return unmapMate;

        boolean unmapSuppAlignment = false;

        if(!unmapRead && checkUnmapSupplementaryAlignments(read))
        {
            clearSupplementaryAlignment(read);
            mStats.SuppAlignmentCount.incrementAndGet();
            unmapSuppAlignment = true;
        }

        return unmapRead || unmapMate || unmapSuppAlignment;
    }

    private boolean checkTransformUnpairedRead(final SAMRecord read, final UnmapRegionState regionState)
    {
        // first check the read's alignment itself
        if(read.getReadUnmappedFlag())
            return false; // nothing to do and there won't be a supplementary

        boolean isSupplementary = read.getSupplementaryAlignmentFlag();

        boolean unmapRead = false;

        UnmapReason unmapReason = checkUnmapRead(read, regionState);

        if(unmapReason != UnmapReason.NONE)
        {
            unmapRead = true;

            switch(unmapReason)
            {
                case HIGH_DEPTH:
                    mStats.HighDepthCount.incrementAndGet();
                    break;

                case SOFT_CLIP:
                    mStats.LongSoftClipCount.incrementAndGet();
                    break;

                default:
                    break;
            }
        }
        else if(isSupplementary)
        {
            // supplementaries are unmapped if their primary, or an associated supplementary, will be unmapped
            unmapRead = checkUnmapSupplementaryRead(read);
        }

        if(unmapRead && isSupplementary)
        {
            // these will be dropped from the BAM
            setUnmappedAttributes(read);
            mStats.SupplementaryCount.incrementAndGet();
            return true;
        }

        if(unmapRead)
        {
            unmapReadAlignment(read, false, true);
            mStats.UnmappedCount.incrementAndGet();
        }

        // nothing more to do for a supplementary whose primary has been unmapped, or if it was unmapped
        if(isSupplementary)
            return false;

        boolean unmapSuppAlignment = false;

        if(!unmapRead && checkUnmapSupplementaryAlignments(read))
        {
            clearSupplementaryAlignment(read);
            mStats.SuppAlignmentCount.incrementAndGet();
            unmapSuppAlignment = true;
        }

        return unmapRead || unmapSuppAlignment;
    }

    private boolean checkSecondaryRead(final SAMRecord read, final UnmapRegionState regionState)
    {
        if(read.getReadUnmappedFlag())
            return false; // nothing to do

        UnmapReason unmapReason = checkUnmapRead(read, regionState);

        if(unmapReason == UnmapReason.NONE)
            return false;

        switch(unmapReason)
        {
            case HIGH_DEPTH:
                mStats.HighDepthCount.incrementAndGet();
                break;

            case SOFT_CLIP:
                mStats.LongSoftClipCount.incrementAndGet();
                break;

            default:
                break;
        }

        // these will be dropped from the BAM - only set the unmapped flag
        read.setReadUnmappedFlag(true);
        mStats.SecondaryCount.incrementAndGet();
        return true;
    }

    private enum UnmapReason
    {
        HIGH_DEPTH,
        SOFT_CLIP,
        CHIMERIC,
        NONE;
    }

    private UnmapReason checkUnmapRead(final SAMRecord read, final UnmapRegionState regionState)
    {
        RegionMatchType matchType = findMaxDepthRegionOverlap(
                read.getAlignmentStart(), read.getAlignmentEnd(), regionState.PartitionRegions, regionState, true);

        return checkUnmapRead(read, matchType);
    }

    private UnmapReason checkUnmapRead(final SAMRecord read, final RegionMatchType matchType)
    {
        if(matchType == RegionMatchType.NONE)
            return UnmapReason.NONE;

        if(matchType == RegionMatchType.HIGH_DEPTH)
            return UnmapReason.HIGH_DEPTH;

        if(exceedsSoftClipLength(read.getCigarString()))
            return UnmapReason.SOFT_CLIP;

        if(read.getReadPairedFlag())
        {
            if(isChimericRead(read, false))
                return UnmapReason.CHIMERIC;
        }

        return UnmapReason.NONE;
    }

    private boolean checkUnmapMate(final SAMRecord read, final RegionMatchType matchType)
    {
        if(matchType == RegionMatchType.NONE)
            return false;

        if(matchType == RegionMatchType.HIGH_DEPTH)
            return true;

        if(isChimericRead(read, true))
            return true;

        if(exceedsSoftClipLength(read.getStringAttribute(MATE_CIGAR_ATTRIBUTE)))
            return true;

        return false;
    }

    private boolean checkUnmapSupplementaryRead(final SAMRecord read)
    {
        // a supplementary should be dropped if its primary satisfies the unmapping criteria
        List<SupplementaryReadData> alignments = SupplementaryReadData.extractAlignments(read);

        if(alignments == null)
            return false;

        for(SupplementaryReadData suppData : alignments)
        {
            RegionMatchType matchType = supplementaryMaxDepthRegionOverlap(suppData);

            if(matchType == RegionMatchType.NONE)
                continue;

            if(matchType == RegionMatchType.HIGH_DEPTH)
                return true;

            if(exceedsSoftClipLength(suppData.Cigar))
                return true;

            if(isSupplementaryChimericRead(read, suppData))
                return true;
        }

        return false;
    }

    private boolean checkUnmapSupplementaryAlignments(final SAMRecord read)
    {
        // checks if any supplementary alignment qualifies for unmapping
        List<SupplementaryReadData> alignments = SupplementaryReadData.extractAlignments(read);

        if(alignments == null)
            return false;

        for(SupplementaryReadData suppData : alignments)
        {
            if(checkUnmapSupplementaryAlignment(suppData))
                return true;
        }

        return false;
    }

    private boolean checkUnmapSupplementaryAlignment(final SupplementaryReadData suppData)
    {
        if(suppData == null)
            return false;

        RegionMatchType matchType = supplementaryMaxDepthRegionOverlap(suppData);

        if(matchType == RegionMatchType.NONE)
            return false;

        if(matchType == RegionMatchType.HIGH_DEPTH)
            return true;

        if(exceedsSoftClipLength(suppData.Cigar))
            return true;

        return false;
    }

    private enum RegionMatchType
    {
        NONE,
        HIGH_DEPTH,
        OTHER;
    }

    private static final int READ_END_APPROX_BUFFER = 100;

    public boolean mateInUnmapRegion(final SAMRecord read)
    {
        return mateMaxDepthRegionOverlap(read, null) != RegionMatchType.NONE;
    }

    public RegionMatchType mateMaxDepthRegionOverlap(final SAMRecord read, @Nullable final UnmapRegionState regionState)
    {
        // first check for a local mate vs the partition's unmapped regions
        boolean checkLocalRegions = false;

        List<HighDepthRegion> mateRegions;

        if(regionState != null && read.getMateReferenceName().equals(regionState.Partition.chromosome())
        && regionState.Partition.containsPosition(read.getMateAlignmentStart()))
        {
            mateRegions = regionState.PartitionRegions;
            checkLocalRegions = true;
        }
        else
        {
            mateRegions = mChrLocationsMap.get(read.getMateReferenceName());

            if(mateRegions == null)
                return RegionMatchType.NONE;

            if(!isWithinRegionRange(read.getMateAlignmentStart(), mateRegions))
                return RegionMatchType.NONE;
        }

        // coordinate check must be identical to how the mate checks itself, which requires knowledge of aligned bases
        int mateEnd;

        if(read.hasAttribute(MATE_CIGAR_ATTRIBUTE))
        {
            mateEnd = getReadEndFromCigarStr(read.getMateAlignmentStart(), read.getStringAttribute(MATE_CIGAR_ATTRIBUTE));
        }
        else
        {
            // TODO: this won't work for split reads
            int approxReadEnd = read.getMateAlignmentStart() + read.getReadBases().length + READ_END_APPROX_BUFFER;
            mateEnd = approxReadEnd;
        }

        return findMaxDepthRegionOverlap(
                read.getMateAlignmentStart(), mateEnd, mateRegions,
                checkLocalRegions ? regionState : null, false);
    }

    private RegionMatchType supplementaryMaxDepthRegionOverlap(final SupplementaryReadData suppReadData)
    {
        final List<HighDepthRegion> suppRegions = mChrLocationsMap.get(suppReadData.Chromosome);

        if(suppRegions == null)
            return RegionMatchType.NONE;

        if(!isWithinRegionRange(suppReadData.Position, suppRegions))
            return RegionMatchType.NONE;

        int readEnd = getReadEndFromCigarStr(suppReadData.Position, suppReadData.Cigar);
        return findMaxDepthRegionOverlap(suppReadData.Position, readEnd, suppRegions, null, false);
    }

    private static boolean isWithinRegionRange(final int readStart, final List<HighDepthRegion> regions)
    {
        int startIndex = binarySearch(readStart, regions);
        for(int i = startIndex; i < regions.size(); ++i)
        {
            HighDepthRegion region = regions.get(i);

            if(positionsOverlap(readStart, readStart + READ_END_APPROX_BUFFER, region.start(), region.end()))
                return true;

            if(region.start() > readStart + READ_END_APPROX_BUFFER)
                break;
        }

        return false;
    }

    private static final int NO_INDEX_MATCH = -1;

    private int checkRegionStateMatch(final int readStart, final int readEnd, final UnmapRegionState regionState)
    {
        if(regionState == null || regionState.LastMatchedRegionIndex == null || regionState.PartitionRegions.isEmpty())
            return NO_INDEX_MATCH;

        if(!positionsWithin(readStart, readEnd, regionState.Partition.start(), regionState.Partition.end()))
            return NO_INDEX_MATCH;

        HighDepthRegion region = regionState.PartitionRegions.get(regionState.LastMatchedRegionIndex);

        if(readStart < region.start())
        {
            if(regionState.LastMatchedRegionIndex == 0)
                return regionState.LastMatchedRegionIndex; // returning the first region is still valid since the read is within the partition

            HighDepthRegion prevRegion = regionState.PartitionRegions.get(regionState.LastMatchedRegionIndex - 1);

            return readStart >= prevRegion.start() ? regionState.LastMatchedRegionIndex - 1 : NO_INDEX_MATCH;
        }

        if(region.containsPosition(readStart))
            return regionState.LastMatchedRegionIndex;

        if(regionState.LastMatchedRegionIndex >= regionState.PartitionRegions.size() - 1)
            return NO_INDEX_MATCH;

        HighDepthRegion nextRegion = regionState.PartitionRegions.get(regionState.LastMatchedRegionIndex + 1);

        return readStart < nextRegion.start() ? regionState.LastMatchedRegionIndex : NO_INDEX_MATCH;
    }

    private RegionMatchType findMaxDepthRegionOverlap(
            final int readStart, final int readEnd, final List<HighDepthRegion> regions,
            @Nullable final UnmapRegionState regionState, boolean updateRegionState)
    {
        if(regions.isEmpty())
            return RegionMatchType.NONE;

        RegionMatchType matchType = RegionMatchType.NONE;

        int startIndex = checkRegionStateMatch(readStart, readEnd, regionState);

        if(startIndex == NO_INDEX_MATCH)
        {
            startIndex = binarySearch(readStart, regions);

            if(updateRegionState)
                regionState.LastMatchedRegionIndex = startIndex;
        }

        // in effect the binary search finds a current overlap or the previous region, so at most 2 regions will be tested
        for(int i = startIndex; i < regions.size(); ++i)
        {
            HighDepthRegion region = regions.get(i);
            if(region.start() > readEnd)
                break;

            if(!overlapsRegion(region, readStart, readEnd))
                continue;

            if(!positionsOverlap(readStart, readEnd, region.start(), region.end()))
                continue;

            if(positionsWithin(readStart, readEnd, region.start(), region.end()))
            {
                return region.maxDepth() >= UNMAP_MIN_HIGH_DEPTH ? RegionMatchType.HIGH_DEPTH : RegionMatchType.OTHER;
            }

            int overlapBases = min(region.end(), readEnd) - max(region.start(), readStart) + 1;
            int readLength = readEnd - readStart + 1;
            int nonOverlappingBases = readLength - overlapBases;

            if(nonOverlappingBases < UNMAP_MAX_NON_OVERLAPPING_BASES)
            {
                if(region.maxDepth() >= UNMAP_MIN_HIGH_DEPTH)
                    return RegionMatchType.HIGH_DEPTH;
                else
                    matchType = RegionMatchType.OTHER; // and continue searching for a better match
            }
        }

        return matchType;
    }

    public static boolean overlapsRegion(final HighDepthRegion region, int readStart, int readEnd)
    {
        if(!positionsOverlap(readStart, readEnd, region.start(), region.end()))
            return false;

        if(positionsWithin(readStart, readEnd, region.start(), region.end()))
            return true;

        int overlapBases = min(region.end(), readEnd) - max(region.start(), readStart) + 1;
        int readLength = readEnd - readStart + 1;
        int nonOverlappingBases = readLength - overlapBases;

        return nonOverlappingBases < UNMAP_MAX_NON_OVERLAPPING_BASES;
    }

    private static boolean isChimericRead(final SAMRecord record, boolean checkForMate)
    {
        boolean readUnmapped = record.getReadUnmappedFlag();
        boolean mateUnmapped = record.getMateUnmappedFlag();

        // or a fragment length outside the expected maximum
        if(abs(record.getInferredInsertSize()) > UNMAP_CHIMERIC_FRAGMENT_LENGTH_MAX)
            return true;

        // an unmapped mate
        if(!checkForMate && mateUnmapped)
            return true;

        if(checkForMate && readUnmapped)
            return true;

        if(record.getReadPairedFlag())
        {
            // inter-chromosomal
            if(!record.getReferenceName().equals(record.getMateReferenceName()))
                return true;

            // inversion
            if(record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag())
                return true;
        }

        return false;
    }

    private static boolean exceedsSoftClipLength(final String cigar)
    {
        return getClipLengthFromCigarStr(cigar) >= UNMAP_MIN_SOFT_CLIP;
    }

    private static boolean isSupplementaryChimericRead(final SAMRecord record, SupplementaryReadData suppReadData)
    {
        if(record.getReadPairedFlag())
        {
            if(record.getMateUnmappedFlag())
                return true;

            // inter-chromosomal
            if(!suppReadData.Chromosome.equals(record.getMateReferenceName()))
                return true;

            // inversion
            boolean primaryReadOnNegativeStrand = suppReadData.Strand == SupplementaryReadData.SUPP_NEG_STRAND;

            if(primaryReadOnNegativeStrand == record.getMateNegativeStrandFlag())
                return true;
        }

        return false;
    }

    private static void setUnmappedAttributes(final SAMRecord read)
    {
        read.setReadUnmappedFlag(true);
        read.setSecondaryAlignment(false);
        read.setMappingQuality(0);

        // store the original mapping
        setUnmapCoordsAttribute(read, read.getReferenceName(), read.getAlignmentStart());
    }

    public static void unmapReadAlignment(final SAMRecord read, final boolean mateUnmapped, final boolean unmapMate)
    {
        setUnmappedAttributes(read);

        read.setProperPairFlag(false);

        // clear insert size
        read.setInferredInsertSize(0);

        // clear reference index, reference name and alignment start
        if(mateUnmapped || unmapMate)
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

        // If mate is unmapped, clear its reference name and alignment start here, since no unmapping logic will be run for the mate
        if(mateUnmapped)
        {
            read.setMateAlignmentStart(0);
            read.setMateReferenceIndex(NO_CHROMOSOME_INDEX);
            read.setMateReferenceName(NO_CHROMOSOME_NAME);
        }

        // clear cigar if present
        read.setCigarString(NO_CIGAR);

        // clear the supplementary data too
        clearSupplementaryAlignment(read);
    }

    private static final String UNMAPP_COORDS_DELIM = ":";

    public static String[] parseUnmappedCoords(final String mateCoordsStr)
    {
        return mateCoordsStr.split(UNMAPP_COORDS_DELIM, 2);
    }

    private static void setUnmapCoordsAttribute(final SAMRecord read, final String chromosome, final int position)
    {
        // some alt contigs have the delimiter in their name, hence the replace
        read.setAttribute(
                UNMAP_ATTRIBUTE,
                chromosome.replaceAll(UNMAPP_COORDS_DELIM, "") + UNMAPP_COORDS_DELIM + position);
    }

    public static void unmapMateAlignment(final SAMRecord read, final boolean readUnmapped, final boolean unmapRead)
    {
        // set flag mate unmapped
        read.setMateUnmappedFlag(true);
        read.setProperPairFlag(false);

        // store the mate's original mapping
        setUnmapCoordsAttribute(read, read.getMateReferenceName(), read.getMateAlignmentStart());

        // clear insert size
        read.setInferredInsertSize(0);

        // set mate reference index, reference name and alignment start
        if(readUnmapped || unmapRead)
        {
            // set both to unknown
            read.setMateAlignmentStart(0);
            read.setMateReferenceIndex(NO_CHROMOSOME_INDEX);
            read.setMateReferenceName(NO_CHROMOSOME_NAME);

            // clear the read's coords too - this is not handled elsewhere
            if(readUnmapped)
            {
                read.setAlignmentStart(0);
                read.setReferenceIndex(NO_CHROMOSOME_INDEX);
                read.setReferenceName(NO_CHROMOSOME_NAME);
            }
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

    public static void clearSupplementaryAlignment(final SAMRecord read)
    {
        // remove supplementary alignment details from attributes
        if(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
        {
            read.setAttribute(SUPPLEMENTARY_ATTRIBUTE, null);
        }
    }

    private static Map<String,List<HighDepthRegion>> loadUnmapRegions(final String filename)
    {
        return UnmappedRegions.loadUnmapRegions(filename);
    }

    @VisibleForTesting
    public static int getClipLengthFromCigarStr(final String cigarStr)
    {
        if(cigarStr == null || cigarStr.isEmpty())
            return 0;

        int softClipCount = 0;
        int elementLength = 0;
        for(int i = 0; i < cigarStr.length(); ++i)
        {
            final char c = cigarStr.charAt(i);
            if(c == 'S' || c == 'H')
            {
                softClipCount += elementLength;
                elementLength = 0;
                continue;
            }

            int digit = c - '0';
            if(digit >= 0 && digit <= 9)
            {
                elementLength = elementLength * 10 + digit;
            }
            else
            {
                elementLength = 0;
            }
        }

        return softClipCount;
    }

    private static int getReadEndFromCigarStr(final int readStart, final String cigarStr)
    {
        int currentPosition = readStart;
        int elementLength = 0;

        for(int i = 0; i < cigarStr.length(); ++i)
        {
            final char c = cigarStr.charAt(i);
            final boolean isAddItem = (c == 'D' || c == 'M' || c == 'N');

            if(isAddItem)
            {
                currentPosition += elementLength;
                elementLength = 0;
                continue;
            }

            final int digit = c - '0';
            if(digit >= 0 && digit <= 9)
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

    @VisibleForTesting
    public void addRegion(final String chromosome, final HighDepthRegion region)
    {
        List<HighDepthRegion> regions = mChrLocationsMap.get(chromosome);
        if(regions == null)
        {
            regions = Lists.newArrayList();
            mChrLocationsMap.put(chromosome, regions);
        }

        regions.add(region);

        mEnabled = true;
    }
}
