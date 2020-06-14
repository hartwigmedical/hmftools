package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.REALIGN_MAX_SOFT_CLIP_BASE_LENGTH;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.REALIGN_MIN_SOFT_CLIP_BASE_LENGTH;

import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.isofox.common.ReadRecord;

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

        /*
        Set<Byte> orientations = Sets.newHashSet();
        reads.stream().filter(x -> !x.isSupplementaryAlignment()).forEach(x -> orientations.add(x.orientation()));
        return orientations.size() == 1;
         */
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


}
