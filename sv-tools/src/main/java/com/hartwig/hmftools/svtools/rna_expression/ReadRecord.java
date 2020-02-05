package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.RegionReadData.extractTransId;

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

    private final Map<RegionReadData,Integer> mMappedRegions; // regions related to this read and their match type

    private final Map<String,Integer> mTranscriptClassification;

    // type of match between a read region (ie a 'M' section in CIGAR) and an exon
    public static final int MATCH_TYPE_NONE = -1;
    public static final int MATCH_TYPE_EXON_BOUNDARY = 0; // read matches one exon boundary
    public static final int MATCH_TYPE_WITHIN_EXON = 1; // read fully contained within the exon
    public static final int MATCH_TYPE_EXON_MATCH = 2; // read fully contained within the exon
    public static final int MATCH_TYPE_UNSPLICED = 3; // reads spanning to unmapped regions where adjacent regions exist
    public static final int MATCH_TYPE_INTRONIC = 4;

    // type of match against a transcript
    public static final int TRANS_MATCH_UNKONWN = -1;
    public static final int TRANS_MATCH_ALT = 0; // maps to unknown exons (ie intronic regions)
    public static final int TRANS_MATCH_UNSPLICED = 1; // has reads spanning exon-intron boundaries
    public static final int TRANS_MATCH_EXONIC = 2; // has reads fully within an exon
    public static final int TRANS_MATCH_SPLICE_JUNCTION = 3; // 2 or more exons correctly covered by a read
    public static final int TRANS_MATCH_OTHER_TRANS = 4; // correctly matched to another trans

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

            int matchType = getRegionMatchType(region);

            mMappedRegions.put(region, matchType);
        }

        for(final String transId : transcripts)
        {
            int transMatchType = TRANS_MATCH_UNKONWN;

            List<RegionReadData> transRegions = regions.stream()
                    .filter(x -> x.getRefRegions().stream().anyMatch(y -> y.contains(transId)))
                    .collect(Collectors.toList());

            // if any reads cross and exon-intron boundary, then mark the transcript as unspliced

            if(transRegions.size() == 1 && mMappedCoords.size() == 1)
            {
                // simple case of a single exon and read section
                RegionReadData region = transRegions.get(0);
                int matchType = mMappedRegions.get(region);

                if (matchType == MATCH_TYPE_NONE)
                {
                    transMatchType = TRANS_MATCH_ALT;
                }
                else if (matchType == MATCH_TYPE_UNSPLICED)
                {
                    transMatchType = TRANS_MATCH_UNSPLICED;
                }
            }
            else if(mMappedCoords.size() > transRegions.size())
            {
                transMatchType = TRANS_MATCH_ALT;
            }
            else
            {
                int minExonRank = 0;
                int maxExonRank = 0;

                for (int regionIndex = 0; regionIndex < transRegions.size(); ++regionIndex)
                {
                    RegionReadData region = transRegions.get(regionIndex);

                    int exonRank = region.getExonRank(transId);
                    maxExonRank = max(maxExonRank, exonRank);
                    minExonRank = minExonRank == 0 ? exonRank : min(exonRank, minExonRank);

                    int mappingIndex = getRegionMappingIndex(region);

                    if (mappingIndex < 0 || mappingIndex != regionIndex)
                    {
                        transMatchType = TRANS_MATCH_ALT;
                        break;
                    }

                    int matchType = mMappedRegions.get(region);

                    if (matchType == MATCH_TYPE_UNSPLICED)
                    {
                        transMatchType = TRANS_MATCH_UNSPLICED;
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
                                transMatchType = TRANS_MATCH_ALT;
                                break;
                            }
                        }
                        else if(regionIndex == transRegions.size() - 1)
                        {
                            if(missStart)
                            {
                                transMatchType = TRANS_MATCH_ALT;
                                break;
                            }
                        }
                        else if(missStart || missEnd)
                        {
                            transMatchType = TRANS_MATCH_ALT;
                            break;
                        }
                    }
                }

                if (transMatchType == TRANS_MATCH_UNKONWN)
                {
                    int expectedRegions = maxExonRank - minExonRank + 1;
                    if (transRegions.size() < expectedRegions)
                        transMatchType = TRANS_MATCH_ALT;
                }
            }

            if(transMatchType == TRANS_MATCH_UNKONWN)
            {
                if(transRegions.size() > 1)
                {
                    transMatchType = TRANS_MATCH_SPLICE_JUNCTION;
                }
                else
                {
                    transMatchType = TRANS_MATCH_EXONIC;
                }
            }

            mTranscriptClassification.put(transId, transMatchType);
        }
    }

    public static boolean validTranscriptType(int transType)
    {
        return transType == TRANS_MATCH_EXONIC || transType == TRANS_MATCH_SPLICE_JUNCTION;
    }

    public static boolean validRegionMatchType(int matchType)
    {
        return matchType == MATCH_TYPE_EXON_BOUNDARY || matchType == MATCH_TYPE_WITHIN_EXON || matchType == MATCH_TYPE_EXON_MATCH;
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

    public int getRegionMatchType(final RegionReadData region)
    {
        int mappingIndex = getRegionMappingIndex(region);
        if (mappingIndex < 0)
            return MATCH_TYPE_NONE;

        return getRegionMatchType(region, mappingIndex);
    }

    private int getRegionMatchType(final RegionReadData region, int mappingIndex)
    {
        if(mappingIndex < 0 || mappingIndex >= mMappedCoords.size())
            return MATCH_TYPE_NONE;

        final long[] readSection = mMappedCoords.get(mappingIndex);
        long readStartPos = readSection[SE_START];
        long readEndPos = readSection[SE_END];

        if (readEndPos < region.start() || readStartPos > region.end())
            return MATCH_TYPE_NONE;

        if (readStartPos < region.start() || readEndPos > region.end())
            return MATCH_TYPE_UNSPLICED;

        if (readStartPos > region.start() && readEndPos < region.end())
            return MATCH_TYPE_WITHIN_EXON;

        if (readStartPos == region.start() && readEndPos == region.end())
            return MATCH_TYPE_EXON_MATCH;

        return MATCH_TYPE_EXON_BOUNDARY;
    }

    // TODO - adjust to set base matching only for matched transcript regions
    public static int setMatchingBases(final RegionReadData region, final ReadRecord read)
    {
        if(read.Cigar == null)
            return MATCH_TYPE_NONE;

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
                long readStartPos = read.PosStart + readBaseIndex - readSkippedBases + regionSkippedBases;

                if(hasSplitRegions)
                {
                    // check whether this matched section corresponds to the region (and not an adjacent exon region)
                    long readPosEnd = readStartPos + element.getLength() - 1;

                    if(readPosEnd < region.start() || readStartPos > region.end())
                    {
                        // this read relates to a difference region
                        readBaseIndex += element.getLength();
                        continue;
                    }
                }

                int matchType = MATCH_TYPE_NONE;

                if(region.start() < readStartPos)
                {
                    regionBaseIndex += readStartPos - region.start();
                }
                else if(region.start() > readStartPos)
                {
                    readBaseIndex += region.start() - readStartPos;
                }

                if(regionBaseIndex < 0 || regionBaseIndex >= regionBases.length())
                {
                    LOGGER.warn("invalid base match: cigar({}) read(pos={} start={} index={}) region(pos={} index={} skipped={}) matchLen({})",
                            cigar.toString(), readStartPos, read.PosStart, readBaseIndex,
                            region.start() + regionBaseIndex, regionBaseIndex, regionSkippedBases, element.getLength());
                    return MATCH_TYPE_NONE;
                }

                // check the bases in this matched region
                int matchedBases = 0;
                int overlapBases = 0;

                for(int i = 0; i < element.getLength(); ++i)
                {
                    if(regionBaseIndex >= regionBases.length() || readBaseIndex >= readBases.length())
                        break;

                    ++overlapBases;

                    if(regionBases.charAt(regionBaseIndex) == readBases.charAt(readBaseIndex))
                    {
                        ++regionMatchedBases[regionBaseIndex];
                        ++matchedBases;
                    }

                    ++regionBaseIndex;
                    ++readBaseIndex;
                }

                // classify the type of match
                long readEndPos = read.PosStart + readBaseIndex - readSkippedBases + regionSkippedBases - 1;

                if(element.getLength() == matchedBases)
                {
                    if (readStartPos == region.start() || readEndPos == region.end())
                    {
                        matchType = MATCH_TYPE_EXON_BOUNDARY;
                    }
                    else
                    {
                        matchType = MATCH_TYPE_WITHIN_EXON;
                    }
                }
                else if(matchedBases == overlapBases)
                {
                    // the read crosses this region boundary but matched for the portion which overlapped
                    matchType = MATCH_TYPE_UNSPLICED;
                }
                else
                {
                    LOGGER.trace("insufficient base match: cigar({}) read(index={} pos={}) region(index={} pos={}) matchLen({}) matched({})",
                            cigar.toString(), readBaseIndex, read.PosStart + readBaseIndex, regionBaseIndex, region.start() + regionBaseIndex,
                            element.getLength(), matchedBases);
                }

                return matchType;
            }
        }

        // not an issue since must map to other regions
        LOGGER.trace("no base match: cigar({}) read({}->{} index={}) region(index={} pos={})",
                cigar.toString(), read.PosStart, read.PosEnd, readBaseIndex,
                region.start(), region.end(), regionBaseIndex);

        return MATCH_TYPE_NONE;
    }


    private static int MIN_BASE_MATCH = 4;

    private boolean checkSoftClippedRegions(final RegionReadData region)
    {
        boolean matched = false;

        // look for soft-clipped bases which match the exon before or afterwards
        if(Cigar.getFirstCigarElement().getOperator() == CigarOperator.S
        && Cigar.getFirstCigarElement().getLength() >= MIN_BASE_MATCH && !region.getPreRegions().isEmpty())
        {
            int scLength = Cigar.getFirstCigarElement().getLength();
            final String scBases = ReadBases.substring(0, scLength);

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

                    //region.addMatchedRead(MATCH_TYPE_SPLICE_JUNCTION);
                    //preRegion.addMatchedRead(MATCH_TYPE_SPLICE_JUNCTION);

                    preRegion.addLinkedRegion(region);
                    region.addLinkedRegion(preRegion);

                    break;
                }
            }
        }

        if(Cigar.getLastCigarElement().getOperator() == CigarOperator.S
        && Cigar.getLastCigarElement().getLength() >= MIN_BASE_MATCH && !region.getPostRegions().isEmpty())
        {
            int scLength = Cigar.getLastCigarElement().getLength();
            final String scBases = ReadBases.substring(ReadBases.length() - scLength, ReadBases.length());

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

                    // region.addMatchedRead(MATCH_TYPE_SPLICE_JUNCTION);
                    // postRegion.addMatchedRead(MATCH_TYPE_SPLICE_JUNCTION);

                    postRegion.addLinkedRegion(region);
                    region.addLinkedRegion(postRegion);

                    break;
                }
            }
        }

        return matched;
    }


    public final Map<RegionReadData,Integer> getMappedRegions() { return mMappedRegions; }

    public final List<RegionReadData> getMappedRegions(int matchType)
    {
        return mMappedRegions.entrySet().stream()
                .filter(x -> x.getValue() == matchType).map(x -> x.getKey()).collect(Collectors.toList());
    }

    public final Map<String,Integer> getTranscriptClassifications() { return mTranscriptClassification; }

    public int getTranscriptClassification(final String trans)
    {
        Integer transType = mTranscriptClassification.get(trans);
        return transType != null ? transType : TRANS_MATCH_UNKONWN;
    }

    public String toString()
    {
        return String.format("range(%s: %d -> %d, range=%d) length(%d) cigar(%s) id(%s)",
                Chromosome, PosStart, PosEnd, range(), Length, Cigar != null ? Cigar.toString() : "", Id);
    }

}
