package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.isofox.common.CommonUtils.canonicalAcceptor;
import static com.hartwig.hmftools.isofox.common.CommonUtils.canonicalDonor;
import static com.hartwig.hmftools.isofox.common.CommonUtils.deriveCommonRegions;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.NONE;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.matchRank;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.REALIGN_MAX_SOFT_CLIP_BASE_LENGTH;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.REALIGN_MIN_SOFT_CLIP_BASE_LENGTH;
import static com.hartwig.hmftools.isofox.fusion.FusionJunctionType.CANONICAL;
import static com.hartwig.hmftools.isofox.fusion.FusionJunctionType.KNOWN;
import static com.hartwig.hmftools.isofox.fusion.FusionJunctionType.UNKNOWN;
import static com.hartwig.hmftools.isofox.fusion.FusionTransExon.fromList;
import static com.hartwig.hmftools.isofox.fusion.FusionTransExon.mergeUnique;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.TransExonRef;

import org.jetbrains.annotations.Nullable;

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
            if(junctOrientations[seIndex] == ORIENT_REV)
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

    public static void checkFusionPositionAdjustmentsVsKnownExons(
            final FusionReadData fusion, final List<TranscriptData>[] transcriptLists, final RefGenomeInterface refGenome)
    {
        int positionBuffer = fusion.splitJunctionOverlap();
        boolean adjustmentsApplied = false;

        int[] requiredPositionAdjusts = new int[] {0, 0};

        for(int se = SE_START; se <= SE_END; ++se)
        {
            List<TranscriptData> transDataList = transcriptLists[se];
            List<TransExonRef> matchedTransExons = Lists.newArrayList();

            ExonBoundaryMatch topExonMatch = null;

            // look for an exact exon boundary match using the overlap and record if it would require a position adjustment
            int juncPosition = fusion.junctionPositions()[se];
            byte juncOrient = fusion.junctionOrientations()[se];

            // take any previous adjustment from the other position's known exon match
            if(requiredPositionAdjusts[se] != 0)
            {
                juncPosition -= requiredPositionAdjusts[se] * juncOrient;
            }

            for(TranscriptData transData : transDataList)
            {
                ExonBoundaryMatch exonMatch = findExonBoundaryMatch(transData, juncPosition, juncOrient, positionBuffer);

                if(exonMatch == null)
                    continue;

                if(topExonMatch == null || (exonMatch.isExact() && !topExonMatch.isExact()))
                {
                    topExonMatch = exonMatch;

                    // only keep transcripts matching the best
                    matchedTransExons.clear();
                    matchedTransExons.add(exonMatch.transExonRef());
                }
                else if(exonMatch.isExact() == topExonMatch.isExact())
                {
                    matchedTransExons.add(exonMatch.transExonRef());
                }
            }

            if(!matchedTransExons.isEmpty())
            {
                // purge other previously added non-matching transcript references
                fusion.getTransExonRefsByPos(se).clear();
                fusion.getTransExonRefsByPos(se).addAll(matchedTransExons);
            }

            if(topExonMatch != null)
            {
                fusion.junctionTypes()[se] = FusionJunctionType.KNOWN;

                if(!adjustmentsApplied)
                {
                    // adjust the positions if required
                    if(topExonMatch.boundaryPosition() != juncPosition)
                    {
                        int positionAdjust = abs(topExonMatch.boundaryPosition() - juncPosition);

                        requiredPositionAdjusts[se] = positionAdjust;

                        // adjust the other position by the remainder
                        int otherSe = switchIndex(se);
                        int otherPositionAdjust = positionBuffer - positionAdjust;
                        requiredPositionAdjusts[otherSe] = otherPositionAdjust;
                    }
                    else
                    {
                        // the other position must absorb any position adjustment
                        int otherSe = switchIndex(se);
                        requiredPositionAdjusts[otherSe] = positionBuffer;
                    }

                    // cancel any check for position adjustment for the upper / other position
                    positionBuffer = 0;
                    adjustmentsApplied = true;
                }
            }
        }

        if(adjustmentsApplied)
        {
            fusion.markSplitJunctionOverlapApplied();

            for(int se = SE_START; se <= SE_END; ++se)
            {
                int positionAdjustment = requiredPositionAdjusts[se];

                if(positionAdjustment != 0)
                {
                    fusion.junctionPositions()[se] -= positionAdjustment * fusion.junctionOrientations()[se];

                    // reset junction bases now the position has shifted
                    fusion.setJunctionBases(refGenome, se);
                }
            }
        }
    }

    public static void checkFusionPositionAdjustmentsVsCanonicalSpliceSites(final FusionReadData fusion, final RefGenomeInterface refGenome)
    {
        // now fusion genes have been identified, check for canonical splice sites, and apply any remaining position adjustments
        checkCanonicalSpliceJunction(fusion);

        if(fusion.splitJunctionOverlap() == 0 || fusion.splitJunctionOverlapApplied())
            return;

        int positionAdjustment = fusion.splitJunctionOverlap();

        // favour known over canonical over unknown
        int[] juncPositions = fusion.junctionPositions();
        byte[] juncOrientations = fusion.junctionOrientations();

        int juncStartTypeOrdinal = fusion.junctionTypes()[SE_START].ordinal();
        int juncEndTypeOrdinal = fusion.junctionTypes()[SE_END].ordinal();

        int[] positionsAdjustments = new int[] {0, 0};

        if(juncStartTypeOrdinal < juncEndTypeOrdinal)
        {
            // for a DEL this shifts the position up further up since orientation is -ve
            positionsAdjustments[SE_END] = positionAdjustment;
        }
        else if(juncEndTypeOrdinal < juncStartTypeOrdinal)
        {
            positionsAdjustments[SE_START] = positionAdjustment;
        }
        else
        {
            // split the change
            int halfOverlap = positionAdjustment / 2;

            positionsAdjustments[SE_START] = (positionAdjustment % 2) == 0 ? halfOverlap : halfOverlap + 1; // round up if an odd length
            positionsAdjustments[SE_END] = positionAdjustment - positionsAdjustments[SE_START];
        }

        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(positionsAdjustments[se] != 0)
            {
                juncPositions[se] -= positionsAdjustments[se] * juncOrientations[se];
                fusion.setJunctionBases(refGenome, se);

                if(fusion.junctionTypes()[se] != KNOWN)
                {
                    // this can revert an existing canonical junction back to unknown after a position adjustment
                    fusion.junctionTypes()[se] = matchesCanonicalSpliceJunction(
                            juncOrientations[se], fusion.junctionSpliceBases()[se], fusion.geneStrandByPosition(se)) ? CANONICAL : UNKNOWN;
                }
            }
        }
    }

    private record ExonBoundaryMatch(TransExonRef transExonRef, int boundaryPosition, boolean isExact) {}

    private static ExonBoundaryMatch findExonBoundaryMatch(
            final TranscriptData transData, final int juncPosition, final byte juncOrientation, final int positionBuffer)
    {
        int juncPosLower = juncPosition;
        int juncPosUpper = juncPosition;

        if(positionBuffer > 0)
        {
            if(juncOrientation == ORIENT_FWD)
                juncPosLower -= positionBuffer;
            else
                juncPosUpper += positionBuffer;
        }

        for(ExonData exon : transData.exons())
        {
            int exonPosition = juncOrientation == ORIENT_FWD ? exon.End : exon.Start;

            if(positionWithin(exonPosition, juncPosLower, juncPosUpper))
            {
                TransExonRef transExonRef = new TransExonRef(
                        transData.GeneId, transData.TransId, transData.TransName, exon.Rank, transData.IsCanonical);

                boolean isExact = juncPosition == exonPosition;

                return new ExonBoundaryMatch(transExonRef, exonPosition, isExact);
            }
        }

        return null;
    }

    public static void checkCanonicalSpliceJunction(final FusionReadData fusion)
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            // now that the strandedness of the fusion has been determined, check for canonical splice sites if not matching known
            if(fusion.junctionTypes()[se] == KNOWN)
                continue;

            if(matchesCanonicalSpliceJunction(
                    fusion.junctionOrientations()[se], fusion.junctionSpliceBases()[se], fusion.geneStrandByPosition(se)))
            {
                fusion.junctionTypes()[se] = CANONICAL;
            }
        }
    }

    public static boolean matchesCanonicalSpliceJunction(
            final byte juncOrientation, final String juncSpliceBases, @Nullable final Byte juncStrand)
    {
        if(juncStrand != null)
        {
            boolean isDonor = (juncOrientation == juncStrand);

            if(isDonor && canonicalDonor(juncSpliceBases, juncStrand))
                return true;
            else if(!isDonor && canonicalAcceptor(juncSpliceBases, juncStrand))
                return true;
        }
        else
        {
            // try them both if strand is not known
            byte asDonorStrand = juncOrientation;
            byte asAcceptorStrand = (byte)(-juncOrientation);

            if(canonicalDonor(juncSpliceBases, asDonorStrand) || canonicalAcceptor(juncSpliceBases, asAcceptorStrand))
                return true;
        }

        return false;
    }
}
