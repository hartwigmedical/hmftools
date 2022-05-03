package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.isofox.common.CommonUtils.deriveCommonRegions;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.NONE;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.matchRank;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.REALIGN_MAX_SOFT_CLIP_BASE_LENGTH;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.REALIGN_MIN_SOFT_CLIP_BASE_LENGTH;
import static com.hartwig.hmftools.isofox.fusion.FusionTransExon.fromList;
import static com.hartwig.hmftools.isofox.fusion.FusionTransExon.mergeUnique;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.TransExonRef;

public class FusionUtils
{
    public static String formLocation(final String chromosome, final int geneCollectionId, boolean isGenic)
    {
        if(isGenic)
            return String.format("%s:%d", chromosome, geneCollectionId);
        else
            return String.format("%s:pre_%d", chromosome, geneCollectionId);
    }

    public static String formChromosomePair(final String[] chr) { return formChromosomePair(chr[SE_START], chr[SE_END]); }
    public static String formChromosomePair(final String chr1, final String chr2) { return chr1 + "_" + chr2; }

    public static int[] findSplitReadJunction(final FusionRead read)
    {
        if(!read.ContainsSplit)
            return null;

        return read.junctionPositions();
    }

    public static boolean hasRealignableSoftClip(final FusionRead read, int se, boolean checkMax)
    {
        return (read.SoftClipLengths[se] >= REALIGN_MIN_SOFT_CLIP_BASE_LENGTH
                && (!checkMax || read.SoftClipLengths[se] <= REALIGN_MAX_SOFT_CLIP_BASE_LENGTH));
    }

    public static boolean isRealignedFragmentCandidate(final FusionRead read)
    {
        return hasRealignableSoftClip(read, SE_START, true) || hasRealignableSoftClip(read, SE_END, true);
    }

    public static void setMaxSplitMappedLength(
            int seIndex, final List<FusionRead> reads, final int[] junctPositions, final byte[] junctOrientations, final int[] maxSplitLengths)
    {
        // find the longest section mapped across the junction
        final List<FusionRead> matchingReads = reads.stream()
                .filter(x -> positionWithin(junctPositions[seIndex], x.posStart(), x.posEnd())).collect(Collectors.toList());

        if(matchingReads.isEmpty()) // can occur with the fragments from a fusion merged in due to homology
            return;

        List<int[]> mappedCoords;

        if(matchingReads.size() == 1)
        {
            mappedCoords = matchingReads.get(0).MappedCoords;
        }
        else
        {
            mappedCoords = deriveCommonRegions(matchingReads.get(0).MappedCoords, matchingReads.get(1).MappedCoords);
        }

        int mappedBases = 0;

        for(int[] coord : mappedCoords)
        {
            if(junctOrientations[seIndex] == NEG_ORIENT)
            {
                if(coord[SE_END] < junctPositions[seIndex])
                    continue;

                mappedBases += coord[SE_END] - max(junctPositions[seIndex], coord[SE_START]) + 1;
            }
            else
            {
                if(coord[SE_START] > junctPositions[seIndex])
                    break;

                mappedBases += min(junctPositions[seIndex], coord[SE_END]) - coord[SE_START] + 1;
            }
        }

        maxSplitLengths[seIndex] = max(mappedBases, maxSplitLengths[seIndex]);
    }

    public static RegionMatchType extractTopTransExonRefs(
            final Map<RegionMatchType,List<TransExonRef>> transExonRefMap,
            final RegionMatchType existingMatchType, final List<FusionTransExon> existingTransExonRefs)
    {
        RegionMatchType topMatchType = NONE;
        List<TransExonRef> topTransExonRefs = null;

        for(Map.Entry<RegionMatchType, List<TransExonRef>> entry : transExonRefMap.entrySet())
        {
            RegionMatchType matchType = entry.getKey();

            if(matchRank(matchType) < matchRank(existingMatchType))
                continue;

            if(topMatchType == null || matchRank(matchType) >= matchRank(topMatchType))
            {
                topMatchType = matchType;
                topTransExonRefs = entry.getValue();
            }
        }

        if(topTransExonRefs == null)
            return existingMatchType;

        if(matchRank(topMatchType) > matchRank(existingMatchType))
            existingTransExonRefs.clear();

        mergeUnique(existingTransExonRefs, fromList(topTransExonRefs));

        return topMatchType;
    }

    public static void checkMissingGeneData(final FusionRead read, final List<TranscriptData> transDataList)
    {
        if(!read.IsGenicRegion[SE_END] || !read.spansGeneCollections())
            return;

        // due to the way the BAM fragment allocator processes reads per gene collection, the upper gene collection will have missed its
        // transcript exon data, so populate this now

        int upperCoordIndex = read.MappedCoords.size() - 1;
        final int[] upperCoords = read.MappedCoords.get(upperCoordIndex);

        List<FusionTransExon> transExonRefs = Lists.newArrayList();
        RegionMatchType topMatchType = NONE;

        for(TranscriptData transData : transDataList)
        {
            if(!positionsWithin(upperCoords[SE_START], upperCoords[SE_END], transData.TransStart, transData.TransEnd))
                continue;

            for(ExonData exonData : transData.exons())
            {
                if(!positionsOverlap(upperCoords[SE_START], upperCoords[SE_END], exonData.Start, exonData.End))
                    continue;

                if(exonData.Start > upperCoords[SE_END])
                    break;

                RegionMatchType matchType;
                if(upperCoords[SE_START] == exonData.Start || upperCoords[SE_END] == exonData.End)
                {
                    matchType = RegionMatchType.EXON_BOUNDARY;
                }
                else if(positionsWithin(upperCoords[SE_START], upperCoords[SE_END], exonData.Start, exonData.End))
                {
                    matchType = RegionMatchType.WITHIN_EXON;
                }
                else
                {
                    matchType = RegionMatchType.EXON_INTRON;
                }

                if(matchRank(matchType) < matchRank(topMatchType))
                    continue;

                if(matchRank(matchType) > matchRank(topMatchType))
                {
                    transExonRefs.clear();
                    topMatchType = matchType;
                }

                FusionTransExon teRef = new FusionTransExon(transData.TransId, exonData.Rank);
                transExonRefs.add(teRef);
            }
        }

        read.setUpperTransExonRefs(transExonRefs, topMatchType);
    }

    public static final String SUPP_ALIGNMENT_DELIM = ",";

    public static Integer suppAlignmentPosition(final String suppAlignment)
    {
        // 21,39794900,-,33M43S,255,0;
        if(suppAlignment == null)
            return null;

        final String[] items = suppAlignment.split(SUPP_ALIGNMENT_DELIM);
        return items.length >= 5 ? Integer.parseInt(items[1]) : null;
    }

    public static String suppAlignmentChromosome(final String suppAlignment)
    {
        if(suppAlignment == null)
            return null;

        final String[] items = suppAlignment.split(SUPP_ALIGNMENT_DELIM);
        return items.length >= 5 ? items[0] : null;
    }

}
