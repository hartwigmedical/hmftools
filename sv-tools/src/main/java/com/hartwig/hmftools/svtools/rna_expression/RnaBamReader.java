package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.MAX_READ_COUNT;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.SAMSlicer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class RnaBamReader
{
    private final RnaExpConfig mConfig;
    private IndexedFastaSequenceFile mIndexedFastaSequenceFile;
    private final File mRefGenomeFile;
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
        mRefGenomeFile = new File(mConfig.RefGenomeFile);

        mReadRecords = Lists.newArrayList();
        mBamRecordCount = 0;
        mCurrentGene = null;

        mSamReader = SamReaderFactory.makeDefault().referenceSequence(mRefGenomeFile).open(new File(mConfig.BamFile));

        try
        {
            LOGGER.debug("loading indexed fasta reference file");
            mIndexedFastaSequenceFile = new IndexedFastaSequenceFile(new File(mConfig.RefGenomeFile));
        }
        catch (IOException e)
        {
            LOGGER.error("Reference file loading failed: {}", e.toString());
            return;
        }
    }

    public void readBamCounts(final GeneReadData geneReadData, final GenomeRegion genomeRegion)
    {
        mReadRecords.clear();
        mBamRecordCount = 0;

        mCurrentGene = geneReadData;

        SAMSlicer samSlicer = new SAMSlicer(DEFAULT_MIN_MAPPING_QUALITY, Lists.newArrayList(genomeRegion));
        samSlicer.slice(mSamReader, this::processSamRecord);
    }

    private void processSamRecord(@NotNull final SAMRecord record)
    {
        // skip records if either end isn't in one of the genic regions
        boolean startMatched = mCurrentGene.getRegionReadData().stream()
                .anyMatch(x -> record.getStart() >= x.Region.start() && record.getStart() <= x.Region.end());

        boolean endMatched = mCurrentGene.getRegionReadData().stream()
                .anyMatch(x -> record.getEnd() >= x.Region.start() && record.getEnd() <= x.Region.end());

        if(!startMatched || !endMatched)
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
            final String regionRefBases = mIndexedFastaSequenceFile.getSubsequenceAt(
                    region.chromosome(), region.start(), region.end()).getBaseString();

            region.setRefBases(regionRefBases);
        }

        // for each record find all exons with an overlap
        // skip records if either end isn't in one of the exons for this gene

        for(final ReadRecord record : mReadRecords)
        {
            // the read covers a part of the exon
            List<RegionReadData> containingRegions = mCurrentGene.getRegionReadData().stream()
                    .filter(x -> x.Region.start() <= record.PosStart && x.Region.end() >= record.PosEnd)
                    .collect(Collectors.toList());

            // the read covers all of the exon (discounting an reads with large insertions or unmapped bases)
            List<RegionReadData> containedRegions = record.range() < 1.1 * record.Length ?
                    mCurrentGene.getRegionReadData().stream()
                    .filter(x -> !containingRegions.contains(x))
                    .filter(x -> x.Region.start() >= record.PosStart && x.Region.end() <= record.PosEnd)
                    .collect(Collectors.toList()) : Lists.newArrayList();

            // the read covers part of the exon and crosses its boundary
            List<RegionReadData> overlappingRegions = mCurrentGene.getRegionReadData().stream()
                    .filter(x -> !containedRegions.contains(x) && !containingRegions.contains(x))
                    .filter(x -> overlaps(x.Region, record))
                    .collect(Collectors.toList());

            if(containedRegions.isEmpty() && containingRegions.isEmpty() && overlappingRegions.isEmpty())
                continue;

            // record the links between these exons
            List<RegionReadData> allRegions = Lists.newArrayList();
            allRegions.addAll(containedRegions);
            allRegions.addAll(containingRegions);
            allRegions.addAll(overlappingRegions);

            if(allRegions.size() > 1)
            {
                for(int i = 0; i < allRegions.size() - 1; ++i)
                {
                    RegionReadData region1 = allRegions.get(i);

                    for(int j = i + 1; j < allRegions.size(); ++j)
                    {
                        RegionReadData region2 = allRegions.get(j);

                        if(region1.getPreRegions().contains(region2) || region1.getPostRegions().contains(region2))
                        {
                            region1.addLinkedRegion(region2);
                            region2.addLinkedRegion(region1);
                        }
                    }
                }
            }
            else
            {
                RegionReadData region = allRegions.get(0);

                // check for a soft-clipping which matches an adjacent exon
                checkSoftClippedRegions(region, record);

                // TODO: also check for exon-intronic sections
            }

            // now check for matching bases in the read vs the reference for the overlapping sections
            allRegions.forEach(x -> setMatchingBases(x, record));
        }
    }

    public static boolean overlaps(final GenomeRegion region, final ReadRecord record)
    {
        // overlapping but neither wholy contained within
        if(region.start() >= record.PosStart && region.start() <= record.PosEnd && region.end() > record.PosEnd)
            return true;

        if(region.end() >= record.PosStart && region.end() <= record.PosEnd && region.start() < record.PosStart)
            return true;

        return false;
    }

    public static void setMatchingBases(final RegionReadData region, final ReadRecord record)
    {
        long overlapStart = max(region.start(), record.PosStart);
        long overlapEnd = min(region.end(), record.PosEnd);
        int overlapLength = (int)(overlapEnd - overlapStart + 1);

        if(overlapLength < 5)
            return;

        // compare the bases at this location
        int recordOverlapStart = (int)(overlapStart - record.PosStart);

        // factor in soft-clipping
        if(record.Cigar.getFirstCigarElement().getOperator() == CigarOperator.S)
        {
            recordOverlapStart += record.Cigar.getFirstCigarElement().getLength();
        }

        int recordOverlapEnd = min(recordOverlapStart + overlapLength, record.ReadBases.length() - 1);

        if(recordOverlapStart < 0 || recordOverlapStart >= record.ReadBases.length())
        {
            recordOverlapStart = 0;
        }

        final String recordBases = record.ReadBases.substring(recordOverlapStart, recordOverlapEnd);

        int regionOverlapStart = (int)(overlapStart - region.start());
        int regionOverlapEnd = regionOverlapStart + overlapLength;
        final String regionBases = region.refBases().substring(regionOverlapStart, regionOverlapEnd);

        if(!recordBases.equals(regionBases))
        {
            // perform a manual comparison
            int matchedBases = findStringOverlaps(regionBases, recordBases);

            if(matchedBases < MIN_BASE_MATCH_PERC * overlapLength)
            {
                LOGGER.trace("region({}) has base-mismatch, overlap({}) matched({}) cigar({})",
                        region, overlapLength, matchedBases, record.Cigar);

                LOGGER.trace("regionBases: pos({} -> {}) {}",
                        regionOverlapStart, regionOverlapEnd, regionBases);
                LOGGER.trace("recordBases: pos({} -> {}) {}",
                        recordOverlapStart, recordOverlapEnd, recordBases);
                return;
            }
        }

        int[] matchedBases = region.refBasesMatched();

        for(long i = overlapStart; i <= overlapEnd; ++i)
        {
            int baseIndex = (int)(i - region.start());
            ++matchedBases[baseIndex];
        }

        region.addMatchedRead();
    }

    private static int MIN_SOFT_CLIPPED_MATCHED_BASES = 5;

    private void checkSoftClippedRegions(final RegionReadData region, final ReadRecord record)
    {
        // look for soft-clipped bases which match the exon before or afterwards
        if(record.Cigar.getFirstCigarElement().getOperator() == CigarOperator.S
        && record.Cigar.getFirstCigarElement().getLength() >= MIN_SOFT_CLIPPED_MATCHED_BASES && !region.getPreRegions().isEmpty())
        {
            int scLength = record.Cigar.getFirstCigarElement().getLength();
            final String scBases = record.ReadBases.substring(0, scLength);
            boolean matched = false;

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

            if(!matched)
                region.addNonAdjacentRead();
        }

        if(record.Cigar.getLastCigarElement().getOperator() == CigarOperator.S
        && record.Cigar.getLastCigarElement().getLength() >= MIN_SOFT_CLIPPED_MATCHED_BASES && !region.getPostRegions().isEmpty())
        {
            int scLength = record.Cigar.getLastCigarElement().getLength();
            final String scBases = record.ReadBases.substring(record.ReadBases.length() - scLength, record.ReadBases.length());
            boolean matched = false;

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

            if(!matched)
                region.addNonAdjacentRead();
        }
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
