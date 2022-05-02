package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.REALIGN_MAX_SOFT_CLIP_BASE_LENGTH;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.REALIGN_MIN_SOFT_CLIP_BASE_LENGTH;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.TransExonRef;

import htsjdk.samtools.CigarOperator;

public final class ChimericUtils
{
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
}
