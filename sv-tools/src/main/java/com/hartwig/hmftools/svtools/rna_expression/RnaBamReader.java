package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.MAX_READ_COUNT;

import java.io.File;
import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.SAMSlicer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class RnaBamReader
{
    private final RnaExpConfig mConfig;
    private final SamReader mSamReader;

    // state relating to the current gene
    private final List<ReadRecord> mReadRecords;
    private int mBamRecordCount;
    private GeneReadData mCurrentGene;

    // private static final int DEFAULT_MIN_BASE_QUALITY = 13;
    private static final int DEFAULT_MIN_MAPPING_QUALITY = 1;

    private static final Logger LOGGER = LogManager.getLogger(RnaBamReader.class);

    public RnaBamReader(final RnaExpConfig config)
    {
        mConfig = config;

        mReadRecords = Lists.newArrayList();
        mBamRecordCount = 0;
        mCurrentGene = null;

        mSamReader = SamReaderFactory.makeDefault().referenceSequence(mConfig.RefGenomeFile).open(new File(mConfig.BamFile));
    }

    public void readBamCounts(final GeneReadData geneReadData, final GenomeRegion genomeRegion)
    {
        mReadRecords.clear();
        mBamRecordCount = 0;

        mCurrentGene = geneReadData;

        SAMSlicer samSlicer = new SAMSlicer(DEFAULT_MIN_MAPPING_QUALITY, Lists.newArrayList(genomeRegion));
        samSlicer.slice(mSamReader, this::processSamRecord);
    }

    public void readBamCounts(final GenomeRegion genomeRegion, final Consumer<SAMRecord> consumer)
    {
        SAMSlicer samSlicer = new SAMSlicer(DEFAULT_MIN_MAPPING_QUALITY, Lists.newArrayList(genomeRegion));
        samSlicer.slice(mSamReader, consumer);
    }

    private void processSamRecord(@NotNull final SAMRecord record)
    {
        // skip records if both ends aren't in one of the genic regions
        boolean startMatched = mCurrentGene.getRegionReadData().stream()
                .anyMatch(x -> record.getStart() >= x.Region.start() && record.getStart() <= x.Region.end());

        boolean endMatched = mCurrentGene.getRegionReadData().stream()
                .anyMatch(x -> record.getEnd() >= x.Region.start() && record.getEnd() <= x.Region.end());

        if(!startMatched && !endMatched)
            return;

        ++mBamRecordCount;

        if(mConfig.ReadCountLimit > 0 && mReadRecords.size() >= mConfig.ReadCountLimit)
            return;

        mReadRecords.add(new ReadRecord(record));
    }

    public void analyseReads()
    {
        mCurrentGene.setTotalReadCount(mBamRecordCount);

        if(mBamRecordCount >= MAX_READ_COUNT)
        {
            LOGGER.warn("gene({}) readCount({}) exceeds max read count", mCurrentGene.GeneData.GeneName, mBamRecordCount);
            // return;
        }

        // cache reference bases for comparison with read bases
        for(RegionReadData region : mCurrentGene.getRegionReadData())
        {
            final String regionRefBases = mConfig.RefFastaSeqFile.getSubsequenceAt(
                    region.chromosome(), region.start(), region.end()).getBaseString();

            region.setRefBases(regionRefBases);
        }

        // for each record find all exons with an overlap
        // skip records if either end isn't in one of the exons for this gene

        for(final ReadRecord read : mReadRecords)
        {
            // the read is fully within the exon
            List<RegionReadData> containingRegions = mCurrentGene.getRegionReadData().stream()
                    .filter(x -> x.Region.start() <= read.PosStart && x.Region.end() >= read.PosEnd)
                    .collect(Collectors.toList());

            // the read covers all of the exon (discounting any intronic or unmapped sections)
            long recordReadLength = read.Length;

            List<RegionReadData> containedRegions = mCurrentGene.getRegionReadData().stream()
                    .filter(x -> !containingRegions.contains(x))
                    .filter(x -> x.length() <= recordReadLength)
                    .filter(x -> x.Region.start() >= read.PosStart && x.Region.end() <= read.PosEnd)
                    .collect(Collectors.toList());

            // the read covers part of the exon and crosses its boundary
            List<RegionReadData> overlappingRegions = mCurrentGene.getRegionReadData().stream()
                    .filter(x -> !containedRegions.contains(x) && !containingRegions.contains(x))
                    .filter(x -> overlaps(x.Region, read))
                    .collect(Collectors.toList());

            if(containedRegions.isEmpty() && containingRegions.isEmpty() && overlappingRegions.isEmpty())
                continue;

            List<RegionReadData> allRegions = Lists.newArrayList();
            allRegions.addAll(containedRegions);
            allRegions.addAll(containingRegions);
            allRegions.addAll(overlappingRegions);

            // now check for matching bases in the read vs the reference for the overlapping sections
            List<RegionReadData> matchedRegions = Lists.newArrayList();
            for(RegionReadData region : allRegions)
            {
                boolean matched = setMatchingBases(region, read);

                if(matched)
                    matchedRegions.add(region);
            }

            boolean linksFound = false;

            if(matchedRegions.size() > 1)
            {
                // record the links between these exons
                List<RegionReadData> linkedRegions = Lists.newArrayList();

                for(int i = 0; i < matchedRegions.size() - 1; ++i)
                {
                    RegionReadData region1 = matchedRegions.get(i);

                    if(linkedRegions.contains(region1))
                        continue;

                    for(int j = i + 1; j < matchedRegions.size(); ++j)
                    {
                        RegionReadData region2 = matchedRegions.get(j);

                        if(region1.getPreRegions().contains(region2) || region1.getPostRegions().contains(region2))
                        {
                            region1.addLinkedRegion(region2);
                            region2.addLinkedRegion(region1);
                            linkedRegions.add(region1);
                            linkedRegions.add(region2);
                        }
                    }
                }

                if(!linkedRegions.isEmpty())
                {
                    linksFound = true;
                    linkedRegions.forEach(x -> matchedRegions.remove(x));
                }
            }

            if(!linksFound && overlappingRegions.size() == 1 && matchedRegions.contains(overlappingRegions.get(0)))
            {
                RegionReadData region = overlappingRegions.get(0);

                // check for an exonic region going into its adjacent intron
                boolean otherRegionMatched = checkNonCodingOverlaps(region, read);

                // check for a soft-clipping which matches an adjacent exon
                if (!otherRegionMatched)
                {
                    checkSoftClippedRegions(region, read);
                }
            }
        }
    }

    public static boolean overlaps(final GenomeRegion region, final ReadRecord record)
    {
        // overlapping but neither wholy contained within
        if(region.start() >= record.PosStart && region.start() <= record.PosEnd && region.end() > record.PosEnd)
            return true; // region starts at or within and ends after

        if(region.end() >= record.PosStart && region.end() <= record.PosEnd && region.start() < record.PosStart)
            return true; // region starts before and ends at the record end or before

        return false;
    }

    public static boolean setMatchingBases(final RegionReadData region, final ReadRecord read)
    {
        if(read.Cigar == null)
        {
            setMatchingBasesNoCigar(region, read);
            return false;
        }

        final Cigar cigar = read.Cigar;

        // first establish whether the read is split across 2 distant regions, and if so which it maps to
        final String readBases = read.ReadBases;
        int readBaseIndex = 0;

        final String regionBases = region.refBases();
        int[] regionMatchedBases = region.refBasesMatched();
        int regionBaseIndex = 0;
        int regionSkippedBases = 0;
        int readSkippedBases = 0;

        boolean hasSplitRegions = cigar.containsOperator(CigarOperator.N)
                && cigar.getCigarElements().stream().filter(x -> x.getOperator() == CigarOperator.M).count() >= 2;

        int matchedRegionsSkipped = 0;

        for(CigarElement element : read.Cigar.getCigarElements())
        {
            if(element.getOperator() == CigarOperator.S || element.getOperator() == CigarOperator.I)
            {
                readBaseIndex += element.getLength();
                readSkippedBases += element.getLength();
            }
            else if(element.getOperator() == CigarOperator.D)
            {
                regionBaseIndex += element.getLength();
            }
            else if(element.getOperator() == CigarOperator.N)
            {
                regionSkippedBases += element.getLength();
            }
            else if(element.getOperator() == CigarOperator.M)
            {
                long readPosition = read.PosStart + readBaseIndex - readSkippedBases + regionSkippedBases;

                if(hasSplitRegions)
                {
                    // check whether this matched section corresponds to the region (and not an adjacent exon region)
                    long readPosEnd = readPosition + element.getLength() - 1;

                    if(readPosEnd < region.start() || readPosition > region.end())
                    {
                        // this read relates to a difference region
                        readBaseIndex += element.getLength();
                        ++matchedRegionsSkipped;
                        continue;
                    }
                }

                if(region.start() < readPosition)
                {
                    regionBaseIndex += readPosition - region.start();
                }
                else if(region.start() > readPosition)
                {
                    readBaseIndex += region.start() - readPosition;
                }

                if(regionBaseIndex < 0 || regionBaseIndex >= regionBases.length())
                {
                    LOGGER.warn("invalid base match: cigar({}) read(pos={} start={} index={}) region(pos={} index={} skipped={}) matchLen({})",
                            cigar.toString(), readPosition, read.PosStart, readBaseIndex,
                            region.start() + regionBaseIndex, regionBaseIndex, regionSkippedBases, element.getLength());
                    return false;
                }

                // check the bases in this matched region
                int matchedBases = 0;

                for(int i = 0; i < element.getLength(); ++i)
                {
                    if(regionBaseIndex >= regionBases.length() || readBaseIndex >= readBases.length())
                        break;

                    if(regionBases.charAt(regionBaseIndex) == readBases.charAt(readBaseIndex))
                    {
                        ++regionMatchedBases[regionBaseIndex];
                        ++matchedBases;
                    }

                    ++regionBaseIndex;
                    ++readBaseIndex;
                }

                if(element.getLength() == matchedBases || matchedBases == region.length()
                || matchedBases >= MIN_BASE_MATCH && matchedBases/(double)element.getLength() > MIN_BASE_MATCH_PERC)
                {
                    region.addMatchedRead();
                }
                else
                {
                    LOGGER.trace("insufficient base match: cigar({}) read(index={} pos={}) region(index={} pos={}) matchLen({}) matched({})",
                            cigar.toString(), readBaseIndex, read.PosStart + readBaseIndex, regionBaseIndex, region.start() + regionBaseIndex,
                            element.getLength(), matchedBases);
                }

                return true;
            }
        }

        // not an issue since must map to other regions
        LOGGER.trace("no base match: cigar({}) read({}->{} index={}) region(index={} pos={})",
                cigar.toString(), read.PosStart, read.PosEnd, readBaseIndex,
                region.start(), region.end(), regionBaseIndex);

        return false;
    }

    private boolean checkNonCodingOverlaps(final RegionReadData region, final ReadRecord read)
    {
        if(read.Cigar == null)
            return false;

        if(read.Cigar.containsOperator(CigarOperator.N))
            return false;

        if (!region.getPreRegions().isEmpty() && read.PosStart < region.start())
        {
            region.addNonAdjacentRead();
            return true;
        }
        else if (!region.getPostRegions().isEmpty() && read.PosEnd > region.end())
        {
            region.addNonAdjacentRead();
            return true;
        }

        return false;
    }

    private static int MIN_BASE_MATCH = 4;

    private boolean checkSoftClippedRegions(final RegionReadData region, final ReadRecord read)
    {
        boolean matched = false;

        // look for soft-clipped bases which match the exon before or afterwards
        if(read.Cigar.getFirstCigarElement().getOperator() == CigarOperator.S
        && read.Cigar.getFirstCigarElement().getLength() >= MIN_BASE_MATCH && !region.getPreRegions().isEmpty())
        {
            int scLength = read.Cigar.getFirstCigarElement().getLength();
            final String scBases = read.ReadBases.substring(0, scLength);

            for(RegionReadData preRegion : region.getPreRegions())
            {
                int baseLength = preRegion.refBases().length();

                if(scBases.length() > baseLength)
                    continue;

                final String endBases = preRegion.refBases().substring(baseLength - scBases.length(), baseLength);
                if(endBases.equals(scBases))
                {
                    matched = true;
                    int[] preRegionBases = preRegion.refBasesMatched();
                    int startBase = baseLength - scBases.length();
                    for(int i = 0; i < scBases.length(); ++i)
                    {
                        ++preRegionBases[startBase + i];
                    }

                    preRegion.addMatchedRead();

                    preRegion.addLinkedRegion(region);
                    region.addLinkedRegion(preRegion);

                    break;
                }
            }
        }

        if(read.Cigar.getLastCigarElement().getOperator() == CigarOperator.S
        && read.Cigar.getLastCigarElement().getLength() >= MIN_BASE_MATCH && !region.getPostRegions().isEmpty())
        {
            int scLength = read.Cigar.getLastCigarElement().getLength();
            final String scBases = read.ReadBases.substring(read.ReadBases.length() - scLength, read.ReadBases.length());

            for(RegionReadData postRegion : region.getPostRegions())
            {
                if(scBases.length() > postRegion.refBases().length())
                    continue;

                final String startBases = postRegion.refBases().substring(0, scBases.length());
                if(startBases.equals(scBases))
                {
                    matched = true;

                    int[] postRegionBases = postRegion.refBasesMatched();
                    for(int i = 0; i < scBases.length(); ++i)
                    {
                        ++postRegionBases[i];
                    }

                    postRegion.addMatchedRead();

                    postRegion.addLinkedRegion(region);
                    region.addLinkedRegion(postRegion);

                    break;
                }
            }
        }

        return matched;
    }

    public static void setMatchingBasesNoCigar(final RegionReadData region, final ReadRecord read)
    {

        // manual base comparison
        // perform a manual comparison
        int matchedBases = findStringOverlaps(region.refBases(), read.ReadBases);

        if(matchedBases < MIN_BASE_MATCH_PERC * region.length())
        {
            /*
            LOGGER.trace("region({}) has base-mismatch, overlap({}) matched({}) cigar({})",
                    region, overlapLength, matchedBases, read.Cigar);

            LOGGER.trace("regionBases: pos({} -> {}) {}",
                    regionBaseStart, regionBaseEnd, regionBases);
            LOGGER.trace("recordBases: pos({} -> {}) {}",
                    recordBaseStart, recordBaseEnd, recordBases);
            return;

         */
        }

        // TODO: support this manual comparison?
        int[] refBasesMatched = region.refBasesMatched();

        /*
        for(long i = overlapStart; i <= overlapEnd; ++i)
        {
            int baseIndex = (int)(i - region.start());
            ++matchedBases[baseIndex];
        }
        */

        region.addMatchedRead();
    }

    private static final double MIN_BASE_MATCH_PERC = 0.9;

    public static int findStringOverlaps(final String str1, final String str2)
    {
        if(str1.length() == 0 || str2.length() == 0)
            return 0;

        int matched = 0;
        int i = 0;
        int j = 0;
        int mismatchIndex = -1;

        // first compare bases at same indices, making note of the first difference if there is one
        while(i < str1.length() && j < str2.length())
        {
            if (str1.charAt(i) == str2.charAt(j))
                ++matched;
            else if(mismatchIndex == -1)
                mismatchIndex = i;

            ++i;
            ++j;
        }

        if(matched > MIN_BASE_MATCH_PERC * min(str1.length(), str2.length()))
            return matched;

        i = j = mismatchIndex;
        matched = mismatchIndex;

        while(i < str1.length() && j < str2.length())
        {
            if(str1.charAt(i) == str2.charAt(j))
            {
                ++i;
                ++j;
                ++matched;
                continue;
            }

            // search ahead in each string in turn for the next short matching sequence
            int startI = i;
            boolean seqFound = false;
            for(; i < str1.length() - 2 && j < str2.length() - 2; ++i)
            {
                if(str1.charAt(i) == str2.charAt(j) && str1.charAt(i+1) == str2.charAt(j+1) && str1.charAt(i+2) == str2.charAt(j+2))
                {
                    seqFound = true;
                    break;
                }
            }

            if(seqFound)
                continue;

            i = startI;

            for(; i < str1.length() - 2 && j < str2.length() - 2; ++j)
            {
                if(str1.charAt(i) == str2.charAt(j) && str1.charAt(i+1) == str2.charAt(j+1) && str1.charAt(i+2) == str2.charAt(j+2))
                {
                    seqFound = true;
                    break;
                }
            }

            if(!seqFound)
                break;
        }

        return matched;
    }

    @VisibleForTesting
    public void addReadRecords(final GeneReadData geneReadData, final List<ReadRecord> readRecords)
    {
        mCurrentGene = geneReadData;
        mBamRecordCount += readRecords.size();
        mReadRecords.clear();;
        mReadRecords.addAll(readRecords);
    }


}
