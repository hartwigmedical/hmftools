package com.hartwig.hmftools.isofox.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_BOUNDARY;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_INTRON;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.WITHIN_EXON;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.exonBoundary;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.validExonMatch;
import static com.hartwig.hmftools.isofox.common.TransMatchType.ALT;
import static com.hartwig.hmftools.isofox.common.TransMatchType.EXONIC;
import static com.hartwig.hmftools.isofox.common.TransMatchType.SPLICE_JUNCTION;
import static com.hartwig.hmftools.isofox.common.TransMatchType.UNKNOWN;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.region.BaseRegion;

public final class ReadTranscriptUtils
{
    private static final int MIN_SC_BASE_MATCH = 2;
    private static final int MAX_SC_BASE_MATCH = 10;
    private static final int MAX_SC_WITHIN_EXON_LENGTH = 2; // must stay below the realignment window (REALIGN_MIN_SOFT_CLIP_BASE_LENGTH)

    public static void processOverlappingRegions(final Read read, final List<RegionReadData> regions)
    {
        MappedCoords mappedCoords = read.mappedCoords();
        Map<RegionReadData,RegionMatchType> mappedRegions = read.getMappedRegions();
        Map<Integer,TransMatchType> transcriptClassification = read.getTranscriptClassifications();

        // process all regions for each transcript as a group to look for inconsistencies with the transcript definition
        Set<Integer> transcripts = Sets.newHashSet();
        boolean hasSoftClipping = read.isLeftClipped() || read.isRightClipped();

        for(RegionReadData region : regions)
        {
            for(TransExonRef ref : region.getTransExonRefs())
            {
                transcripts.add(ref.TransId);
            }

            RegionMatchType matchType = setRegionMatchType(mappedCoords, mappedRegions, region);
            mappedRegions.put(region, matchType);

            boolean checkMissedJunctions = matchType == EXON_INTRON || (hasSoftClipping && exonBoundary(matchType));

            if(checkMissedJunctions)
                checkMissedJunctions(read, mappedCoords, mappedRegions, region);
        }

        int mappedRegionCount = read.mappedRegionCount();

        for(int transId : transcripts)
        {
            // determine for each transcript whether the mapped regions support a spliced transcript, unspliced or alternate splicing
            TransMatchType transMatchType = UNKNOWN;

            List<RegionReadData> transRegions = regions.stream()
                    .filter(x -> x.getTransExonRefs().stream().anyMatch(y -> y.TransId == transId))
                    .collect(Collectors.toList());

            // if any reads cross and exon-intron boundary, then mark the transcript as unspliced

            if(transRegions.size() == 1 && mappedRegionCount == 1)
            {
                // simple case of a single exon and read section
                RegionReadData region = transRegions.get(0);
                RegionMatchType matchType = mappedRegions.get(region);

                if(matchType == RegionMatchType.NONE)
                {
                    // should never happen since implies this read didn't hit the region at all
                    transMatchType = ALT;
                }
                else if(matchType == RegionMatchType.EXON_INTRON)
                {
                    transMatchType = TransMatchType.UNSPLICED;
                }
            }
            else if(mappedRegionCount > transRegions.size())
            {
                transMatchType = ALT;
            }
            else
            {
                int minExonRank = 0;
                int maxExonRank = 0;

                Collections.sort(transRegions);

                for(int regionIndex = 0; regionIndex < transRegions.size(); ++regionIndex)
                {
                    RegionReadData region = transRegions.get(regionIndex);

                    int exonRank = region.getExonRank(transId);
                    maxExonRank = max(maxExonRank, exonRank);
                    minExonRank = minExonRank == 0 ? exonRank : min(exonRank, minExonRank);

                    int mappingIndex = mappedCoords.findRegionIndex(region);
                    int adjustedMappingIndex = mappingIndex;

                    if(mappedCoords.inferredAlignmentAdded(SE_START))
                        --adjustedMappingIndex;

                    if(adjustedMappingIndex < 0 || adjustedMappingIndex != regionIndex)
                    {
                        transMatchType = ALT;
                        break;
                    }

                    RegionMatchType matchType = mappedRegions.get(region);

                    if(matchType == RegionMatchType.EXON_INTRON)
                    {
                        transMatchType = TransMatchType.UNSPLICED;
                        break;
                    }
                    else
                    {
                        BaseRegion readSection = mappedCoords.regionByIndex(mappingIndex);
                        int readStartPos = readSection.start();
                        int readEndPos = readSection.end();

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

                if(transMatchType == UNKNOWN)
                {
                    int expectedRegions = maxExonRank - minExonRank + 1;
                    if(transRegions.size() < expectedRegions)
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

            // any read with soft-clipping which cannot be mapped to the next exon, other than for short likely adapter sequence reads,
            // is classified as alt, unless the clip is within an exon and a threshold
            if(validTranscriptType(transMatchType) && read.containsSoftClipping() && !read.likelyAdaperSoftClipping())
            {
                int softClipSide = read.isLeftClipped() ? SE_START : SE_END;

                if(!mappedCoords.isSoftClipRegionMatched(softClipSide) && !shortClipWithinExon(read, softClipSide, transRegions))
                {
                    transMatchType = ALT;
                }
            }

            transcriptClassification.put(transId, transMatchType);
        }
    }

    private static void checkMissedJunctions(
            final Read read, final MappedCoords mappedCoords, final Map<RegionReadData,RegionMatchType> mappedRegions,
            final RegionReadData region)
    {
        if(read.hasSuppAlignment())
            return;

        // check for reads either soft-clipped or seemingly unspliced, where the extra bases can match with the next exon

        // check start of read
        BaseRegion readSection = mappedCoords.lowestAlignment(false);
        int readStartPos = readSection.start();
        int readEndPos = readSection.end();

        int extraBaseLength = 0;
        int scLength = 0;

        boolean hasRegionOverhang = region.start() > readStartPos && readEndPos > region.start()
                && region.start() - readStartPos <= MAX_SC_BASE_MATCH;

        if(hasRegionOverhang)
        {
            extraBaseLength = region.start() - readStartPos;
        }

        if(read.isLeftClipped() && readStartPos <= region.start())
        {
            scLength = read.leftClipLength();
            extraBaseLength += scLength;
        }

        // less any deleted bases
        // extraBaseLength = max(extraBaseLength - deletedLength, 0);

        // allow a single base match if only 1 region matches
        if(extraBaseLength >= 1 && extraBaseLength <= MAX_SC_BASE_MATCH && scLength <= MAX_SC_BASE_MATCH)
        {
            // first check for a match with the next exon on the lower side
            String extraBases = read.readBases().substring(0, extraBaseLength);

            List<RegionReadData> matchedRegions = region.getPreRegions().stream()
                    .filter(x -> matchesOtherRegionBases(extraBases, x, false)).collect(Collectors.toList());

            if(!matchedRegions.isEmpty())
            {
                mappedCoords.addSoftClipRegionMatched(true, matchedRegions.size());
                mappedRegions.put(region, EXON_BOUNDARY);

                if(matchedRegions.size() == 1 || (matchedRegions.size() > 1 && extraBaseLength < MIN_SC_BASE_MATCH))
                {
                    // truncate the read positions back to match the exon boundary
                    if(!mappedCoords.lowerInferredAlignmentAdded() && hasRegionOverhang)
                        readSection.setStart(readSection.start() + region.start() - readStartPos);
                }

                // if only one region is matched or the min bases matched is satisfied, then create a mapping to the next region,
                // otherwise treat the splice support as ambiguous (it not mapped to the next region)
                if(matchedRegions.size() == 1 || (matchedRegions.size() > 1 && extraBaseLength >= MIN_SC_BASE_MATCH))
                {
                    for(RegionReadData preRegion : matchedRegions)
                    {
                        // add matched coordinates for this exon and add it as a region
                        mappedRegions.put(preRegion, EXON_BOUNDARY);
                        mappedCoords.addInferredRegion(true, preRegion.end() - extraBaseLength + 1, preRegion.end());
                    }
                }
            }
        }

        // check end of read
        readSection = mappedCoords.highestAlignment(false);
        readStartPos = readSection.start();
        readEndPos = readSection.end();

        extraBaseLength = 0;
        scLength = 0;

        hasRegionOverhang = readEndPos > region.end() && readStartPos < region.end() && readEndPos - region.end() <= MAX_SC_BASE_MATCH;

        if(hasRegionOverhang)
        {
            extraBaseLength = readEndPos - region.end();
        }

        if(read.isRightClipped() && readEndPos >= region.end())
        {
            scLength = read.rightClipLength();
            extraBaseLength += scLength;
        }

        if(extraBaseLength >= 1 && extraBaseLength <= MAX_SC_BASE_MATCH && scLength <= MAX_SC_BASE_MATCH)
        {
            // now check for a match to the next exon up
            int readLength = read.baseLength();
            String extraBases = read.readBases().substring(readLength - extraBaseLength, readLength);

            List<RegionReadData> matchedRegions = region.getPostRegions().stream()
                    .filter(x -> matchesOtherRegionBases(extraBases, x, true)).collect(Collectors.toList());

            if(!matchedRegions.isEmpty())
            {
                mappedCoords.addSoftClipRegionMatched(false, matchedRegions.size());

                mappedRegions.put(region, EXON_BOUNDARY);

                if(matchedRegions.size() == 1 || (matchedRegions.size() > 1 && extraBaseLength < MIN_SC_BASE_MATCH))
                {
                    if(!mappedCoords.upperInferredAlignmentAdded() && hasRegionOverhang)
                        readSection.setEnd(readSection.end() - (readEndPos - region.end()));
                }

                if(matchedRegions.size() == 1 || (matchedRegions.size() > 1 && extraBaseLength >= MIN_SC_BASE_MATCH))
                {
                    for(RegionReadData postRegion : matchedRegions)
                    {
                        mappedRegions.put(postRegion, EXON_BOUNDARY);
                        mappedCoords.addInferredRegion(false, postRegion.start(), postRegion.start() + extraBaseLength - 1);
                    }
                }
            }
        }
    }

    private static boolean matchesOtherRegionBases(final String extraBases, final RegionReadData otherRegion, boolean matchToStart)
    {
        int otherRegionLength = otherRegion.length();

        if(extraBases.length() > otherRegionLength)
            return false;

        String otherRegionBases = matchToStart ? otherRegion.refBases().substring(0, extraBases.length())
                : otherRegion.refBases().substring(otherRegionLength - extraBases.length(), otherRegionLength);

        return (otherRegionBases.equals(extraBases));
    }

    private static boolean shortClipWithinExon(final Read read, int se, final List<RegionReadData> transRegions)
    {
        int clipLength = se == SE_START ? read.leftClipLength() : read.rightClipLength();

        if(clipLength > MAX_SC_WITHIN_EXON_LENGTH)
            return false;

        int readPositionBoundary = read.getCoordsBoundary(se);

        int transcriptBoundary = se == SE_START ?
                transRegions.stream().mapToInt(RegionReadData::start).min().orElse(0) :
                transRegions.stream().mapToInt(RegionReadData::end).max().orElse(0);

        return se == SE_START ?
                readPositionBoundary - clipLength >= transcriptBoundary : readPositionBoundary + clipLength <= transcriptBoundary;
    }

    public static List<RegionReadData> getUniqueValidRegion(final Read read1, final Read read2)
    {
        List<RegionReadData> regions = read1.getMappedRegions().entrySet().stream()
                .filter(x -> validExonMatch(x.getValue()))
                .map(x -> x.getKey()).collect(Collectors.toList());

        List<RegionReadData> regions2 = read2.getMappedRegions().entrySet().stream()
                .filter(x -> validExonMatch(x.getValue()))
                .map(x -> x.getKey()).collect(Collectors.toList());

        for(RegionReadData region : regions2)
        {
            if(!regions.contains(region))
                regions.add(region);
        }

        return regions;
    }

    public static boolean validTranscriptType(TransMatchType transType)
    {
        return transType == EXONIC || transType == SPLICE_JUNCTION;
    }

    private static RegionMatchType setRegionMatchType(
            final MappedCoords mappedCoords, final Map<RegionReadData,RegionMatchType> mappedRegions, final RegionReadData region)
    {
        int mappingIndex = mappedCoords.findRegionIndex(region);
        if(mappingIndex == MappedCoords.INVALID_INDEX)
            return RegionMatchType.NONE;

        RegionMatchType matchType = getRegionMatchType(mappedCoords, region, mappingIndex);
        mappedRegions.put(region, matchType);
        return matchType;
    }

    /*
    public static RegionMatchType getRegionMatchType(final MappedCoords mappedCoords, final RegionReadData region)
    {
        int mappingIndex = mappedCoords.findRegionIndex(region);
        if(mappingIndex == MappedCoords.INVALID_INDEX)
            return RegionMatchType.NONE;

        return getRegionMatchType(mappedCoords, region, mappingIndex);
    }
    */

    public static RegionMatchType getRegionMatchType(final MappedCoords mappedCoords, final RegionReadData region, int mappingIndex)
    {
        if(mappingIndex == MappedCoords.INVALID_INDEX || mappingIndex >= mappedCoords.alignmentCount())
            return RegionMatchType.NONE;

        BaseRegion readSection = mappedCoords.regionByIndex(mappingIndex);
        int readStartPos = readSection.start();
        int readEndPos = readSection.end();

        if(readEndPos < region.start() || readStartPos > region.end())
            return RegionMatchType.NONE;

        if(readStartPos < region.start() || readEndPos > region.end())
            return RegionMatchType.EXON_INTRON;

        if(readStartPos > region.start() && readEndPos < region.end())
            return WITHIN_EXON;

        return EXON_BOUNDARY;
    }

    public static void markRegionBases(final List<BaseRegion> readCoords, final RegionReadData region)
    {
        int[] regionBaseDepth = region.refBasesMatched();

        if(regionBaseDepth == null)
            return;

        for(BaseRegion readSection : readCoords)
        {
            int readStartPos = readSection.start();
            int readEndPos = readSection.end();

            if(readStartPos > region.end() || readEndPos < region.start())
                continue;

            // process this overlap
            int regionBaseIndex = readStartPos > region.start() ? readStartPos - region.start() : 0;
            int overlap = min(readEndPos, region.end()) - max(readStartPos, region.start()) + 1;

            if(regionBaseIndex + overlap > regionBaseDepth.length)
            {
                ISF_LOGGER.error("region({}) read coords({} -> {}) regionBaseIndex({}) overlap({}) regionLength({})",
                        region, readStartPos, readEndPos, regionBaseIndex, overlap, regionBaseDepth.length);
                return;
            }

            for(int j = regionBaseIndex; j < regionBaseIndex + overlap; ++j)
            {
                ++regionBaseDepth[j];
            }
        }
    }

    public static int calcFragmentLength(final TranscriptData transData, final Read read1, final Read read2)
    {
        int minReadPos = min(read1.alignmentStart(), read2.alignmentStart());
        int maxReadPos = max(read1.alignmentEnd(), read2.alignmentEnd());
        return calcFragmentLength(transData, minReadPos, maxReadPos);
    }

    public static int calcFragmentLength(final TranscriptData transData, final int minReadPos, final int maxReadPos)
    {
        // calculate fragment length within this transcript assuming it has been spliced
        int transcriptBases = 0;
        boolean startFound = false;

        for(ExonData exon : transData.exons())
        {
            if(!startFound)
            {
                if(minReadPos < exon.Start - MAX_SC_BASE_MATCH)
                    break;

                if(minReadPos > exon.End)
                    continue;

                if(maxReadPos <= exon.End)
                {
                    // within same exon
                    return maxReadPos - minReadPos + 1;
                }

                startFound = true;
                transcriptBases = exon.End - max(exon.Start, minReadPos) + 1;
            }
            else
            {
                if(maxReadPos > exon.End)
                {
                    transcriptBases += exon.baseLength();
                }
                else if(maxReadPos < exon.Start)
                {
                    break;
                }
                else
                {
                    transcriptBases += maxReadPos - exon.Start + 1;
                    break;
                }
            }
        }

        return transcriptBases;
    }
}
