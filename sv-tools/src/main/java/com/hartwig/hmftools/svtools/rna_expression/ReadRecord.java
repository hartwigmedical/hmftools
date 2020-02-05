package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.RegionMatchType.EXON_BOUNDARY;
import static com.hartwig.hmftools.svtools.rna_expression.RegionMatchType.EXON_INTRON;
import static com.hartwig.hmftools.svtools.rna_expression.RegionMatchType.EXON_MATCH;
import static com.hartwig.hmftools.svtools.rna_expression.RegionMatchType.WITHIN_EXON;
import static com.hartwig.hmftools.svtools.rna_expression.RegionReadData.extractTransId;
import static com.hartwig.hmftools.svtools.rna_expression.TransMatchType.ALT;
import static com.hartwig.hmftools.svtools.rna_expression.TransMatchType.EXONIC;
import static com.hartwig.hmftools.svtools.rna_expression.TransMatchType.SPLICE_JUNCTION;
import static com.hartwig.hmftools.svtools.rna_expression.TransMatchType.UNKNOWN;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class ReadRecord
{
    private final SAMRecord mSamRecord;

    public final String Id;
    public final String Chromosome;
    public final long PosStart;
    public final long PosEnd;

    public final String ReadBases;
    public final int Length; // of bases
    public final Cigar Cigar;

    public final List<long[]> mMappedCoords;
    private boolean mLowerInferredAdded;
    private boolean mUpperInferredAdded;

    private final Map<RegionReadData,RegionMatchType> mMappedRegions; // regions related to this read and their match type

    private final Map<String,TransMatchType> mTranscriptClassification;

    private static final Logger LOGGER = LogManager.getLogger(ReadRecord.class);

    public static ReadRecord from(final SAMRecord record)
    {
        return new ReadRecord(
                record, record.getReadName(), record.getReferenceName(), record.getStart(), record.getEnd(),
                record.getReadString(), record.getCigar());
    }

    public ReadRecord(
            final SAMRecord record, final String id, final String chromosome, long posStart, long posEnd,
            final String readBases, @NotNull final Cigar cigar)
    {
        mSamRecord = record;

        Id = id;
        Chromosome = chromosome;
        PosStart = posStart;
        PosEnd = posEnd;
        ReadBases = readBases;
        Length = ReadBases.length();
        Cigar = cigar;

        mMappedCoords = Lists.newArrayList();
        generateMappedCoords();

        mMappedRegions = Maps.newHashMap();
        mTranscriptClassification = Maps.newHashMap();
        mLowerInferredAdded = false;
        mUpperInferredAdded = false;
    }

    public final SAMRecord samRecord() { return mSamRecord; }

    public long range() { return PosEnd - PosStart; }

    public List<long[]> getMappedRegionCoords() { return mMappedCoords; }

    public boolean overlapsMappedReads(long posStart, long posEnd)
    {
        return mMappedCoords.stream().anyMatch(x -> (!(x[SE_START] > posEnd || x[SE_END] < posStart)));
    }

    private void generateMappedCoords()
    {
        // first establish whether the read is split across 2 distant regions, and if so which it maps to
        // int readBaseIndex = 0;
        // int readSkippedBases = 0;
        int posOffset = 0;

        for(CigarElement element : Cigar.getCigarElements())
        {
            if(element.getOperator() == CigarOperator.S || element.getOperator() == CigarOperator.I)
            {
                // readBaseIndex += element.getLength();
                // readSkippedBases += element.getLength();
            }
            else if(element.getOperator() == CigarOperator.D)
            {
                // regionBaseIndex += element.getLength();
            }
            else if(element.getOperator() == CigarOperator.N)
            {
                // regionSkippedBases += element.getLength();
                posOffset += element.getLength();
            }
            else if(element.getOperator() == CigarOperator.M)
            {
                // long readStartPos = PosStart + readBaseIndex - readSkippedBases + regionSkippedBases;
                // long readEndPos = PosStart + readBaseIndex - readSkippedBases + regionSkippedBases - 1;
                long readStartPos = PosStart + posOffset;
                long readEndPos = readStartPos + element.getLength() - 1;

                mMappedCoords.add(new long[] {readStartPos, readEndPos});

                posOffset += element.getLength();
            }
        }
    }

    public void processOverlappingRegions(final List<RegionReadData> regions)
    {
        // process all regions for each transcript as a group to look for inconsistencies with the transcript definition
        List<String> transcripts = Lists.newArrayList();

        for(RegionReadData region : regions)
        {
            for(final String ref : region.getRefRegions())
            {
                final String transId = extractTransId(ref);

                if (!transcripts.contains(transId))
                    transcripts.add(transId);
            }

            RegionMatchType matchType = getRegionMatchType(region);
            mMappedRegions.put(region, matchType);

            if(matchType == EXON_INTRON || Cigar.containsOperator(CigarOperator.S))
                checkMissedJunctions(region);
        }

        for(final String transId : transcripts)
        {
            TransMatchType transMatchType = UNKNOWN;

            List<RegionReadData> transRegions = regions.stream()
                    .filter(x -> x.getRefRegions().stream().anyMatch(y -> y.contains(transId)))
                    .collect(Collectors.toList());

            // if any reads cross and exon-intron boundary, then mark the transcript as unspliced

            if(transRegions.size() == 1 && mMappedCoords.size() == 1)
            {
                // simple case of a single exon and read section
                RegionReadData region = transRegions.get(0);
                RegionMatchType matchType = mMappedRegions.get(region);

                if (matchType == RegionMatchType.NONE)
                {
                    transMatchType = ALT;
                }
                else if (matchType == RegionMatchType.EXON_INTRON)
                {
                    transMatchType = TransMatchType.UNSPLICED;
                }
            }
            else if(mMappedCoords.size() > transRegions.size())
            {
                transMatchType = ALT;
            }
            else
            {
                int minExonRank = 0;
                int maxExonRank = 0;

                Collections.sort(transRegions);

                for (int regionIndex = 0; regionIndex < transRegions.size(); ++regionIndex)
                {
                    RegionReadData region = transRegions.get(regionIndex);

                    int exonRank = region.getExonRank(transId);
                    maxExonRank = max(maxExonRank, exonRank);
                    minExonRank = minExonRank == 0 ? exonRank : min(exonRank, minExonRank);

                    int mappingIndex = getRegionMappingIndex(region);

                    if (mappingIndex < 0 || mappingIndex != regionIndex)
                    {
                        transMatchType = ALT;
                        break;
                    }

                    RegionMatchType matchType = mMappedRegions.get(region);

                    if (matchType == RegionMatchType.EXON_INTRON)
                    {
                        transMatchType = TransMatchType.UNSPLICED;
                        break;
                    }
                    else
                    {
                        final long[] readSection = mMappedCoords.get(mappingIndex);
                        long readStartPos = readSection[SE_START];
                        long readEndPos = readSection[SE_END];

                        boolean missStart = readStartPos > region.start();
                        boolean missEnd = readEndPos < region.end();
                        if(regionIndex == 0)
                        {
                            if(missEnd)
                            {
                                transMatchType = ALT;
                                break;
                            }
                        }
                        else if(regionIndex == transRegions.size() - 1)
                        {
                            if(missStart)
                            {
                                transMatchType = ALT;
                                break;
                            }
                        }
                        else if(missStart || missEnd)
                        {
                            transMatchType = ALT;
                            break;
                        }
                    }
                }

                if (transMatchType == UNKNOWN)
                {
                    int expectedRegions = maxExonRank - minExonRank + 1;
                    if (transRegions.size() < expectedRegions)
                        transMatchType = ALT;
                }
            }

            if(transMatchType == UNKNOWN)
            {
                if(transRegions.size() > 1)
                {
                    transMatchType = SPLICE_JUNCTION;
                }
                else
                {
                    transMatchType = EXONIC;
                }
            }

            mTranscriptClassification.put(transId, transMatchType);
        }
    }

    public static boolean validTranscriptType(TransMatchType transType)
    {
        return transType == EXONIC || transType == SPLICE_JUNCTION;
    }

    public static boolean validRegionMatchType(RegionMatchType matchType)
    {
        return matchType == EXON_BOUNDARY || matchType == WITHIN_EXON || matchType == EXON_MATCH;
    }

    public int getRegionMappingIndex(final RegionReadData region)
    {
        for(int i = 0; i < mMappedCoords.size(); ++i)
        {
            final long[] readSection = mMappedCoords.get(i);
            long readStartPos = readSection[SE_START];
            long readEndPos = readSection[SE_END];

            if (!(readStartPos > region.end() || readEndPos < region.start()))
                return i;
        }

        return -1;
    }

    public RegionMatchType getRegionMatchType(final RegionReadData region)
    {
        int mappingIndex = getRegionMappingIndex(region);
        if (mappingIndex < 0)
            return RegionMatchType.NONE;

        return getRegionMatchType(region, mappingIndex);
    }

    private RegionMatchType getRegionMatchType(final RegionReadData region, int mappingIndex)
    {
        if(mappingIndex < 0 || mappingIndex >= mMappedCoords.size())
            return RegionMatchType.NONE;

        final long[] readSection = mMappedCoords.get(mappingIndex);
        long readStartPos = readSection[SE_START];
        long readEndPos = readSection[SE_END];

        if (readEndPos < region.start() || readStartPos > region.end())
            return RegionMatchType.NONE;

        if (readStartPos < region.start() || readEndPos > region.end())
            return RegionMatchType.EXON_INTRON;

        if (readStartPos > region.start() && readEndPos < region.end())
            return WITHIN_EXON;

        if (readStartPos == region.start() && readEndPos == region.end())
            return EXON_MATCH;

        return EXON_BOUNDARY;
    }

    public void markRegionBases(final RegionReadData region)
    {
        for(int i = 0; i < mMappedCoords.size(); ++i)
        {
            final long[] readSection = mMappedCoords.get(i);
            long readStartPos = readSection[SE_START];
            long readEndPos = readSection[SE_END];

            if (readStartPos > region.end() || readEndPos < region.start())
                continue;

            // process this overlap
            int regionBaseIndex = readStartPos > region.start() ? (int)(readStartPos - region.start()) : 0;
            int overlap = (int)(min(readEndPos, region.end()) - max(readStartPos, region.start())) + 1;

            int[] regionBaseDepth = region.refBasesMatched();

            if(regionBaseIndex + overlap >= regionBaseDepth.length)
            {
                LOGGER.error("read({}) region({}) coords({} -> {}) regionBaseIndex({}) overlap({})",
                        this.toString(), region, readStartPos, readEndPos, regionBaseIndex, overlap);
                return;
            }

            for(int j = regionBaseIndex; j <= overlap; ++j)
            {
                ++regionBaseDepth[j];
            }

            return;
        }
    }

    private static int MIN_BASE_MATCH = 4;

    private void checkMissedJunctions(final RegionReadData region)
    {
        // check for reads either soft-clipped or apparently unspliced, where the extra bases can match with the next exon

        // check start of read
        long[] readSection = mLowerInferredAdded ? mMappedCoords.get(1) : mMappedCoords.get(0);
        long readStartPos = readSection[SE_START];
        long readEndPos = readSection[SE_END];

        int extraBaseLength = 0;

        if(region.start() > readStartPos && readEndPos > region.start())
        {
            extraBaseLength = (int)(region.start() - readStartPos);
        }

        if(Cigar.getFirstCigarElement().getOperator() == CigarOperator.S)
        {
            extraBaseLength += Cigar.getFirstCigarElement().getLength();
        }

        if(extraBaseLength >= MIN_BASE_MATCH)
        {
            final String extraBases = ReadBases.substring(0, extraBaseLength);

            for(RegionReadData preRegion : region.getPreRegions())
            {
                int baseLength = preRegion.length();

                if(extraBases.length() > baseLength)
                    continue;

                final String endBases = preRegion.refBases().substring(baseLength - extraBases.length(), baseLength);
                if(endBases.equals(extraBases))
                {
                    // add matched coordinates for this exon and it as a region
                    mMappedRegions.put(region, EXON_BOUNDARY);
                    mMappedRegions.put(preRegion, EXON_BOUNDARY);
                    addInferredMappingRegion(true, preRegion.end() - extraBaseLength + 1, preRegion.end());
                }
            }
        }

        readSection = mUpperInferredAdded ? mMappedCoords.get(mMappedCoords.size() - 2) : mMappedCoords.get(mMappedCoords.size() - 1);
        readStartPos = readSection[SE_START];
        readEndPos = readSection[SE_END];

        extraBaseLength = 0;

        if(readEndPos > region.end() && readStartPos < region.end())
        {
            extraBaseLength = (int)(readEndPos - region.end());
        }

        if(Cigar.getLastCigarElement().getOperator() == CigarOperator.S)
        {
            extraBaseLength += Cigar.getLastCigarElement().getLength();
        }

        if(extraBaseLength >= MIN_BASE_MATCH)
        {
            final String extraBases = ReadBases.substring(ReadBases.length() - extraBaseLength, ReadBases.length());

            for(RegionReadData postRegion : region.getPostRegions())
            {
                int baseLength = postRegion.length();

                if(extraBases.length() > baseLength)
                    continue;

                final String endBases = postRegion.refBases().substring(0, extraBases.length());
                if(endBases.equals(extraBases))
                {
                    // add matched coordinates for this exon and it as a region
                    mMappedRegions.put(region, EXON_BOUNDARY);
                    mMappedRegions.put(postRegion, EXON_BOUNDARY);

                    addInferredMappingRegion(false, postRegion.start(), postRegion.start() + extraBaseLength - 1);
                }
            }
        }
    }


    private void addInferredMappingRegion(boolean isLower, long posStart, long posEnd)
    {
        if(isLower)
        {
            if (!mLowerInferredAdded)
            {
                mLowerInferredAdded = true;
                mMappedCoords.add(0, new long[] { posStart, posEnd });
            }
            else
            {
                // lengthen the new region if required
                long[] newSection = mMappedCoords.get(0);
                newSection[SE_START] = min(newSection[SE_START], posStart);
            }
        }
        else
        {
            if(!mUpperInferredAdded)
            {
                mUpperInferredAdded = true;
                mMappedCoords.add(new long[] {posStart, posEnd});
            }
            else
            {
                long[] newSection = mMappedCoords.get(mMappedCoords.size() - 1);
                newSection[SE_END] = max(newSection[SE_END], posEnd);
            }
        }
    }

    public final Map<RegionReadData,RegionMatchType> getMappedRegions() { return mMappedRegions; }

    public final List<RegionReadData> getMappedRegions(RegionMatchType matchType)
    {
        return mMappedRegions.entrySet().stream()
                .filter(x -> x.getValue() == matchType).map(x -> x.getKey()).collect(Collectors.toList());
    }

    public final Map<String,TransMatchType> getTranscriptClassifications() { return mTranscriptClassification; }

    public TransMatchType getTranscriptClassification(final String trans)
    {
        TransMatchType transType = mTranscriptClassification.get(trans);
        return transType != null ? transType : UNKNOWN;
    }

    public String toString()
    {
        return String.format("range(%s: %d -> %d, range=%d) length(%d) cigar(%s) id(%s)",
                Chromosome, PosStart, PosEnd, range(), Length, Cigar != null ? Cigar.toString() : "", Id);
    }

}
