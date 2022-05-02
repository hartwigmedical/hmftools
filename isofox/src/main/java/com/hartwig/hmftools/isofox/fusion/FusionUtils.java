package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.isofox.common.CommonUtils.deriveCommonRegions;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.REALIGN_MAX_SOFT_CLIP_BASE_LENGTH;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.REALIGN_MIN_SOFT_CLIP_BASE_LENGTH;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.TransExonRef;

import htsjdk.samtools.CigarOperator;

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

    public static boolean isInversion(final List<ReadRecord> reads)
    {
        // an inversion must a) be same chromosome b) have supplementary alignment c) have same orientations around the chimeric junction
        if(!reads.stream().anyMatch(x -> x.hasSuppAlignment()) || reads.size() != 3)
            return false;

        byte existingOrient = 0;
        String existingChromosome = "";

        for(ReadRecord read : reads)
        {
            if (!read.hasSuppAlignment())
                continue;

            if(existingChromosome.equals(""))
                existingChromosome = read.Chromosome;
            else if(!existingChromosome.equals(read.Chromosome))
                return false;

            int scLeft = read.isSoftClipped(SE_START) ? read.Cigar.getFirstCigarElement().getLength() : 0;
            int scRight = read.isSoftClipped(SE_END) ? read.Cigar.getLastCigarElement().getLength() : 0;

            if (scLeft == 0 && scRight == 0)
                return false;

            byte orientation = scLeft >= scRight ? POS_ORIENT : NEG_ORIENT;

            if(existingOrient == 0)
                existingOrient = orientation;
            else if(existingOrient == orientation)
                return true;
        }

        return false;
    }

    public static int[] findSplitReadJunction(final ReadRecord read)
    {
        if(!read.containsSplit())
            return null;

        final int maxSplitLength = read.Cigar.getCigarElements().stream()
                .filter(x -> x.getOperator() == CigarOperator.N)
                .mapToInt(x -> x.getLength()).max().orElse(0);

        final List<int[]> mappedCoords = read.getMappedRegionCoords();
        for(int i = 0; i < mappedCoords.size() - 1; ++i)
        {
            final int[] lowerCoords = mappedCoords.get(i);
            final int[] upperCoords = mappedCoords.get(i + 1);

            if(upperCoords[SE_START] - lowerCoords[SE_END] - 1 == maxSplitLength)
            {
                return new int[] { lowerCoords[SE_END], upperCoords[SE_START] };
            }
        }

        return null;
    }

    public static int[] findSplitReadJunction(final FusionRead read)
    {
        if(!read.ContainsSplit)
            return null;

        return read.junctionPositions();
    }

    public static boolean hasRealignableSoftClip(final ReadRecord read, int se, boolean checkMax)
    {
        if(!read.isSoftClipped(se))
            return false;

        int scLength = se == SE_START ? read.Cigar.getFirstCigarElement().getLength() : read.Cigar.getLastCigarElement().getLength();

        return (scLength >= REALIGN_MIN_SOFT_CLIP_BASE_LENGTH && (!checkMax || scLength <= REALIGN_MAX_SOFT_CLIP_BASE_LENGTH));
    }

    public static boolean isRealignedFragmentCandidate(final ReadRecord read)
    {
        return hasRealignableSoftClip(read, SE_START, true) || hasRealignableSoftClip(read, SE_END, true);
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

    public static boolean setHasMultipleKnownSpliceGenes(final List<ReadRecord> reads, final List<String[]> knownPairGeneIds)
    {
        // determines whether reads definitely span multiple genes based on the presence of known splice sites at the fragment junction
        ReadRecord splitRead = null;

        final List<TransExonRef>[] junctionTransRefs = new List[] { Lists.newArrayList(), Lists.newArrayList()};

        for(ReadRecord read : reads)
        {
            if(read.containsSplit())
            {
                // find the largest N-split to mark the junction
                final int[] splitJunction = findSplitReadJunction(read);

                if(splitJunction != null)
                {
                    junctionTransRefs[SE_START].addAll(read.getJunctionMatchingTransRefs(splitJunction[SE_START], true));
                    junctionTransRefs[SE_END].addAll(read.getJunctionMatchingTransRefs(splitJunction[SE_END], false));
                    splitRead = read;
                }

                break;
            }

            if(hasRealignableSoftClip(read, SE_START, false))
                junctionTransRefs[SE_START].addAll(read.getJunctionMatchingTransRefs(read.getCoordsBoundary(SE_START), false));

            if(hasRealignableSoftClip(read, SE_END, false))
                junctionTransRefs[SE_END].addAll(read.getJunctionMatchingTransRefs(read.getCoordsBoundary(SE_END), true));
        }

        if(junctionTransRefs[SE_START].isEmpty() || junctionTransRefs[SE_END].isEmpty())
            return false;

        // if this junction supports any single transcript then it's not a chimeric candidate
        boolean matchesKnownPair = false;
        boolean hasGeneMatch = false;

        for(TransExonRef transExonRefStart : junctionTransRefs[SE_START])
        {
            for(TransExonRef transExonRefEnd : junctionTransRefs[SE_END])
            {
                if(transExonRefStart.TransId == transExonRefEnd.TransId)
                {
                    hasGeneMatch = true;

                    if(transExonRefStart.ExonRank == transExonRefEnd.ExonRank - 1) // supports a transcript's known splice site
                        return false;

                    continue;
                }

                if(transExonRefStart.GeneId.equals(transExonRefEnd.GeneId))
                {
                    hasGeneMatch = true;
                }
                else if(!matchesKnownPair)
                {
                    // look for a known fusion pair
                    if(knownPairGeneIds.stream().anyMatch(x -> x[FS_UP].equals(transExonRefStart.GeneId)
                            && x[FS_DOWN].equals(transExonRefEnd.GeneId)))
                    {
                        matchesKnownPair = true;
                    }
                }
            }
        }

        if(!hasGeneMatch || matchesKnownPair)
        {
            if(splitRead != null)
                splitRead.setHasInterGeneSplit(); // will make use of this when handling fusions

            return true;
        }

        return false;
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

    public static void checkMissingGeneData(final FusionRead read, final List<TranscriptData> transDataList)
    {
        if(!read.IsGenicRegion[SE_END])
            return;

        // due to the way the BAM fragment allocator processes reads per gene collection, the upper gene collection will have missed its
        // transcript exon data, so populate this now

        int upperCoordIndex = read.MappedCoords.size() - 1;
        final int[] upperCoords = read.MappedCoords.get(upperCoordIndex);
        final Map<RegionMatchType,List<TransExonRef>> transExonRefMap = read.getTransExonRefs(SE_END);

        for(TranscriptData transData : transDataList)
        {
            if(!positionsWithin(upperCoords[SE_START], upperCoords[SE_END], transData.TransStart, transData.TransEnd))
                continue;

            for(ExonData exonData : transData.exons())
            {
                if(!positionsOverlap(upperCoords[SE_START], upperCoords[SE_END], exonData.Start, exonData.End))
                    continue;

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

                TransExonRef teRef = new TransExonRef(transData.GeneId, transData.TransId, transData.TransName, exonData.Rank);

                final List<TransExonRef> transExonRefs = transExonRefMap.get(matchType);

                if(transExonRefs == null)
                    transExonRefMap.put(matchType, Lists.newArrayList(teRef));
                else
                    transExonRefs.add(teRef);

                break;
            }
        }
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
