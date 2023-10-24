package com.hartwig.hmftools.markdups.common;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CHROMOSOME_INDEX;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.UNMAP_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.mateNegativeStrand;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.mateUnmapped;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.utils.file.FileDelimiters;

import htsjdk.samtools.SAMRecord;

// TODO: Note that MarkDups does not currently consider contigs outside the main Human Chromosomes, e.g. decoys. Note that currently we only unmap
// regions on the main Human Chromosomes. However, this may lead to inconsistencies when running this on non-prod reference genomes,
// e.g. hg38 + decoys. For example, if a primary read is mapped to Chr1 and its mate is mapped to a decoy, and we are unmapping the primary
// read, then the mate coords of the mate on the decoy will not be correctly updated.
public class ReadUnmapper
{
    private final Map<String, List<HighDepthRegion>> mChrLocationsMap; // keyed by chromosome start
    private boolean mEnabled;
    private final UnmapStats mStats;

    @VisibleForTesting
    public static final int MAX_NON_OVERLAPPING_BASES = 9;
    private static final int MIN_SOFT_CLIP = 11;
    @VisibleForTesting
    public static final int MIN_MAX_DEPTH = 1001;
    private static final int CHIMERIC_FRAGMENT_LENGTH_MAX = 1000;

    public ReadUnmapper(final String filename)
    {
        this(loadUnmapRegions(filename));

        if(mEnabled)
        {
            MD_LOGGER.info("loaded {} unmapping regions from {}",
                    mChrLocationsMap.values().stream().mapToInt(x -> x.size()).sum(), filename);
        }
    }

    public ReadUnmapper(final Map<String, List<HighDepthRegion>> chrLocationsMap)
    {
        mChrLocationsMap = chrLocationsMap;
        mEnabled = mChrLocationsMap != null && !mChrLocationsMap.isEmpty();
        mStats = new UnmapStats();
    }

    public List<HighDepthRegion> getRegions(final String chromosome)
    {
        return mChrLocationsMap.get(chromosome);
    }

    public boolean enabled() { return mEnabled; }

    public UnmapStats stats() { return mStats; }

    public boolean checkTransformRead(final SAMRecord read, final List<HighDepthRegion> readRegions)
    {
        if(!mEnabled || readRegions == null)
            return false;

        // first check the read's alignment itself
        final boolean readUnmapped = read.getReadUnmappedFlag();
        final boolean mateUnmapped = mateUnmapped(read);

        if(readUnmapped && mateUnmapped)
            return false; // nothing to do and there won't be a supplementary

        final boolean isSupplementary = read.getSupplementaryAlignmentFlag();

        boolean unmapRead = false;
        if(!readUnmapped)
        {
            if(isSupplementary)
            {
                // Note that we unmap the supplementary if the primary read (contained in the supplementary read's supplementary attribute)
                // is unmapped.
                unmapRead = checkIfUnmapSupplementary(read, true) || checkIfUnmapRead(read, readRegions);
            }
            else
            {
                unmapRead = checkIfUnmapRead(read, readRegions);
            }
        }

        if(unmapRead && isSupplementary)
        {
            // these will be dropped from the BAM
            setUnmappedAttributes(read);
            mStats.SupplementaryCount.incrementAndGet();
            return true;
        }

        final boolean unmapMate = !mateUnmapped && checkIfUnmapMate(read);

        if(unmapRead)
        {
            unmapReadAlignment(read, mateUnmapped, unmapMate);
        }

        if(unmapMate)
        {
            // Note that for supplementary reads, at this point readUnmapped will be false, since a supplementary must be mapped. Also,
            // unmapRead is false, else we would have returned already. Therefore, this will only modify the mate properties and attributes.
            // TODO: This will set the mate coords to this supplementary's coords, not the primary's coords, is this correct?
            unmapMateAlignment(read, readUnmapped, unmapRead);
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

        // Check supplementary.
        // Note that if this is a supplementary read, and its supplementary (i.e. its associated primary read) is being dropped, we will
        // have returned early above.
        if(isSupplementary)
        {
            return unmapMate;
        }

        boolean unmapSuppAlignment = false;
        if(!unmapRead && checkIfUnmapSupplementary(read, false))
        {
            unmapSupplementary(read);
            mStats.SuppAlignmentCount.incrementAndGet();
            unmapSuppAlignment = true;
        }

        return unmapRead || unmapMate || unmapSuppAlignment;
    }

    private static boolean checkIfUnmapRead(final SAMRecord read, final List<HighDepthRegion> readRegions)
    {
        final int maxDepthOverlap = readMaxDepthRegionOverlap(read.getAlignmentStart(), read.getAlignmentEnd(), readRegions);
        if(maxDepthOverlap < 0)
        {
            return false;
        }

        if(maxDepthOverlap >= MIN_MAX_DEPTH)
        {
            return true;
        }

        if(getSoftClipCountFromCigarStr(read.getCigarString()) >= MIN_SOFT_CLIP)
        {
            return true;
        }

        if(isChimericRead(read, false))
        {
            return true;
        }

        return false;
    }

    private boolean checkIfUnmapMate(final SAMRecord read)
    {
        final int maxDepthOverlap = mateMaxDepthRegionOverlap(read);
        if(maxDepthOverlap < 0)
        {
            return false;
        }

        if(maxDepthOverlap >= MIN_MAX_DEPTH)
        {
            return true;
        }

        // Can only check soft clips if the mate cigar attribute is set.
        if(read.hasAttribute(MATE_CIGAR_ATTRIBUTE))
        {
            final String mateCigar = read.getStringAttribute(MATE_CIGAR_ATTRIBUTE);
            if(getSoftClipCountFromCigarStr(mateCigar) >= MIN_SOFT_CLIP)
            {
                return true;
            }
        }

        if(isChimericRead(read, true))
        {
            return true;
        }

        return false;
    }

    private boolean checkIfUnmapSupplementary(final SAMRecord read, final boolean representsPrimaryRead)
    {
        final SupplementaryReadData suppData = SupplementaryReadData.firstAlignmentFrom(read);
        if(suppData == null)
        {
            return false;
        }

        final int readLength = read.getReadBases().length;
        final int maxDepthOverlap = supplementaryMaxDepthRegionOverlap(suppData, readLength);
        if(maxDepthOverlap < 0)
        {
            return false;
        }

        if(maxDepthOverlap >= MIN_MAX_DEPTH)
        {
            return true;
        }

        if(getSoftClipCountFromCigarStr(suppData.Cigar) >= MIN_SOFT_CLIP)
        {
            return true;
        }

        // We don't want to drop supplementaries based on chimeric read criteria, since this is a criteria for a primary read and its
        // primary mate. However, when checking if we are unmapping a supplementary based on whether its primary is unmapped, we do check
        // the chimeric read criteria.
        if(representsPrimaryRead && isSupplementaryChimericRead(read))
        {
            return true;
        }

        return false;
    }

    @VisibleForTesting
    public int mateMaxDepthRegionOverlap(final SAMRecord read)
    {
        // Returns -1 if there is no overlap.
        final List<HighDepthRegion> mateRegions = mChrLocationsMap.get(read.getMateReferenceName());

        if(mateRegions == null)
        {
            return -1;
        }

        // Coordinate check must be identical to how the mate checks itself, which requires knowledge of aligned bases.
        if(read.hasAttribute(MATE_CIGAR_ATTRIBUTE))
        {
            final int approxMateEnd = read.getMateAlignmentStart() + read.getReadBases().length - 1;
            return readMaxDepthRegionOverlap(
                    read.getMateAlignmentStart(), approxMateEnd, read.getStringAttribute(MATE_CIGAR_ATTRIBUTE), mateRegions);
        }
        else
        {
            // TODO: this won't work for split reads
            final int readBaseLength = read.getReadBases().length;
            return readMaxDepthRegionOverlap(
                    read.getMateAlignmentStart(), read.getMateAlignmentStart() + readBaseLength - 1, mateRegions);
        }
    }

    private int supplementaryMaxDepthRegionOverlap(final SupplementaryReadData suppReadData, final int readLength)
    {
        // Returns -1 if there is no overlap.
        final List<HighDepthRegion> suppRegions = mChrLocationsMap.get(suppReadData.Chromosome);
        if(suppRegions == null)
        {
            return -1;
        }

        final int approxSuppEnd = suppReadData.Position + readLength - 1;
        return readMaxDepthRegionOverlap(suppReadData.Position, approxSuppEnd, suppReadData.Cigar, suppRegions);
    }

    private static int readMaxDepthRegionOverlap(final int readStart, final int readEnd, final List<HighDepthRegion> regions)
    {
        // Returns -1 if there is no overlap.
        int maxDepth = -1;
        for(final HighDepthRegion region : regions)
        {
            if(!positionsOverlap(readStart, readEnd, region.start(), region.end()))
                continue;

            if(positionsWithin(readStart, readEnd, region.start(), region.end()))
            {
                maxDepth = max(maxDepth, region.maxDepth());
                continue;
            }

            final int overlapBases = min(region.end(), readEnd) - max(region.start(), readStart) + 1;
            final int readLength = readEnd - readStart + 1;
            final int nonOverlappingBases = readLength - overlapBases;
            if(nonOverlappingBases <= MAX_NON_OVERLAPPING_BASES)
            {
                maxDepth = max(maxDepth, region.maxDepth());
            }
        }

        return maxDepth;
    }

    private static int readMaxDepthRegionOverlap(final int readStart, int approxReadEnd, final String cigar,
            final List<HighDepthRegion> regions)
    {
        // Returns -1 if there is no overlap.

        // Mate and supplementary coords don't have alignment end so this needs to be determined from the cigar, but only
        // do this is if there is a possible overlap with a region.
        int maxDepth = -1;
        int readEnd = -1;
        for(final HighDepthRegion region : regions)
        {
            if(!positionsOverlap(readStart, approxReadEnd, region.start(), region.end()))
                continue;

            // now test with precise read end alignment
            if(readEnd < 0)
            {
                readEnd = getReadEndFromCigarStr(readStart, cigar);
            }

            final int overlapBases = min(region.end(), readEnd) - max(region.start(), readStart) + 1;
            final int readLength = readEnd - readStart + 1;
            final int nonOverlappingBases = readLength - overlapBases;
            if(nonOverlappingBases <= MAX_NON_OVERLAPPING_BASES)
            {
                maxDepth = max(maxDepth, region.maxDepth());
            }
        }

        return maxDepth;
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
    public static int getSoftClipCountFromCigarStr(final String cigarStr)
    {
        int softClipCount = 0;
        int elementLength = 0;
        for(int i = 0; i < cigarStr.length(); ++i)
        {
            final char c = cigarStr.charAt(i);
            if(c == 'S')
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

    private static boolean isChimericRead(final SAMRecord record, boolean checkForMate)
    {
        final boolean readUnmapped = record.getReadUnmappedFlag();
        final boolean mateUnmapped = mateUnmapped(record);

        if(!checkForMate && readUnmapped)
        {
            return false;
        }

        if(checkForMate && mateUnmapped)
        {
            return false;
        }

        // or a fragment length outside the observed distribution
        if(abs(record.getInferredInsertSize()) > CHIMERIC_FRAGMENT_LENGTH_MAX)
        {
            return true;
        }

        // an unmapped mate
        if(!checkForMate && mateUnmapped)
        {
            return true;
        }

        if(checkForMate && readUnmapped)
        {
            return true;
        }

        if(record.getReadPairedFlag())
        {
            // inter-chromosomal
            if(!record.getReferenceName().equals(record.getMateReferenceName()))
            {
                return true;
            }

            // inversion
            if(record.getReadNegativeStrandFlag() == mateNegativeStrand(record))
            {
                return true;
            }
        }

        return false;
    }

    private static boolean isSupplementaryChimericRead(final SAMRecord record)
    {
        // This looks for whether the supplementary's associated primary, whose information is contained in the supplementary read's
        // supplementary attribute, form a chimeric read with the associate primary mate.
        // Note that the associated primary of a supplementary is always initially mapped.
        final SupplementaryReadData suppReadData = SupplementaryReadData.firstAlignmentFrom(record);

        // TODO: Don't check the inferred insert size. Not much we can do here, since the inferred insert size of a supplementary is always
        // based on the supplemartary's coords, and not the associated primary.

        // an unmapped mate
        if(mateUnmapped(record))
        {
            return true;
        }

        if(record.getReadPairedFlag())
        {
            // inter-chromosomal
            final String primaryChromosome = suppReadData.Chromosome;
            if(!primaryChromosome.equals(record.getMateReferenceName()))
            {
                return true;
            }

            // inversion
            final boolean primaryReadOnNegativeStrand = suppReadData.Strand == SupplementaryReadData.SUPP_NEG_STRAND;
            if(primaryReadOnNegativeStrand == mateNegativeStrand(record))
            {
                return true;
            }
        }

        return false;
    }

    private static void setUnmappedAttributes(final SAMRecord read)
    {
        read.setReadUnmappedFlag(true);
        read.setMappingQuality(0);

        // store the original mapping
        setUnmapCoordsAttribute(read, read.getReferenceName(), read.getAlignmentStart());
    }

    public static void unmapReadAlignment(final SAMRecord read, final boolean mateUnmapped, final boolean unmapRead)
    {
        setUnmappedAttributes(read);

        read.setProperPairFlag(false);

        // clear insert size
        read.setInferredInsertSize(0);

        // clear reference index, reference name and alignment start
        if(mateUnmapped || unmapRead)
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

        // If mate is unmapped, have to clear its reference name and alignment start here, since no unmapping logic will be run for the
        // mate.
        if(mateUnmapped)
        {
            read.setMateAlignmentStart(0);
            read.setMateReferenceIndex(NO_CHROMOSOME_INDEX);
            read.setMateReferenceName(NO_CHROMOSOME_NAME);
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
        }
        else
        {
            // set to primary's details
            read.setMateAlignmentStart(read.getAlignmentStart());
            read.setMateReferenceIndex(read.getReferenceIndex());
            read.setMateReferenceName(read.getReferenceName());
        }

        // If read is unmapped, have to clear its reference name and alignment start here, since no unmapping logic will be run for the
        // read.
        if(readUnmapped)
        {
            read.setAlignmentStart(0);
            read.setReferenceIndex(NO_CHROMOSOME_INDEX);
            read.setReferenceName(NO_CHROMOSOME_NAME);
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

    private static Map<String, List<HighDepthRegion>> loadUnmapRegions(final String filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));
            final String delim = FileDelimiters.inferFileDelimiter(filename);

            // Remove header.
            lines.remove(0);

            final Map<String, List<HighDepthRegion>> chrLocationsMap = Maps.newHashMap();
            for(String line : lines)
            {
                final String[] values = line.split(delim, -1);

                final String chromosome = values[0];
                List<HighDepthRegion> regions = chrLocationsMap.get(chromosome);
                if(regions == null)
                {
                    regions = Lists.newArrayList();
                    chrLocationsMap.put(chromosome, regions);
                }

                final int posStart = Integer.parseInt(values[1]);
                final int posEnd = Integer.parseInt(values[2]);
                final int maxDepth = Integer.parseInt(values[3]);

                final HighDepthRegion region = new HighDepthRegion(posStart, posEnd, maxDepth);
                regions.add(region);
            }

            return chrLocationsMap;
        }
        catch(IOException e)
        {
            MD_LOGGER.error("failed to read high-depth regions file {}: {}", filename, e.toString());
            System.exit(1);
            // Unreachable.
            return null;
        }
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
