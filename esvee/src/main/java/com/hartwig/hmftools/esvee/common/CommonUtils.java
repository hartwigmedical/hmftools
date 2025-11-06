package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_N_BYTE;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_INDEL_MAX_GAP;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_INDEL_MAX_OVERLAP;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.redux.BaseQualAdjustment;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

public final class CommonUtils
{
    public static boolean aboveMinQual(byte qual) { return BaseQualAdjustment.aboveLowBaseQual(qual); }
    public static boolean belowMinQual(byte qual) { return BaseQualAdjustment.isLowBaseQual(qual); }

    public static boolean isHighBaseQual(byte qual) { return BaseQualAdjustment.isHighBaseQual(qual, SvConstants.SEQUENCING_TYPE); }
    public static boolean isMediumBaseQual(byte qual) { return BaseQualAdjustment.isMediumBaseQual(qual, SvConstants.SEQUENCING_TYPE); }

    public static boolean isHigherBaseQualCategory(byte qual1, byte qual2)
    {
        if(qual1 == qual2)
            return false;

        boolean lowBq1 = BaseQualAdjustment.isLowBaseQual(qual1);
        boolean lowBq2 = BaseQualAdjustment.isLowBaseQual(qual2);

        if(lowBq1 != lowBq2)
            return lowBq2;

        boolean highBq1 = isHighBaseQual(qual1);
        boolean highBq2 = isHighBaseQual(qual2);

        if(highBq1 != highBq2)
            return highBq1;

        // both medium
        return qual1 > qual2;
    }

    public static boolean isDiscordantFragment(
            final SAMRecord read, final int maxConcordantFragmentLength, @Nullable final SupplementaryReadData suppData)
    {
        if(read.getReadUnmappedFlag() || !read.getReadPairedFlag() || read.getMateUnmappedFlag())
            return false;

        // supplementaries need to check their primary read chromosomes, not their own
        if(read.getSupplementaryAlignmentFlag() && suppData != null)
        {
            if(!suppData.Chromosome.equals(read.getMateReferenceName()))
                return true;
        }
        else if(!read.getReferenceName().equals(read.getMateReferenceName()))
        {
            return true;
        }

        // inversion
        if(read.getReadNegativeStrandFlag() == read.getMateNegativeStrandFlag())
            return true;

        // check fragment length vs the observed max concordant length
        int fragmentSize = abs(read.getInferredInsertSize());

        if(fragmentSize == 0 || (maxConcordantFragmentLength > 0 && fragmentSize >= maxConcordantFragmentLength))
            return true;

        // lastly look for duplication orientation fragments which aren't overlapping fragments
        if(isDuplicationFragment(read, fragmentSize))
            return true;

        return false;
    }

    public static boolean isDuplicationFragment(final SAMRecord read, int fragmentSize)
    {
        if(fragmentSize > 2 * read.getReadBases().length)
        {
            if(!read.getReadNegativeStrandFlag())
            {
                // expect the +ve orientation read to have a lower position
                if(read.getMateAlignmentStart() < read.getAlignmentStart())
                    return true;
            }
            else
            {
                if(read.getAlignmentStart() < read.getMateAlignmentStart())
                    return true;
            }
        }

        return false;
    }

    public static boolean isLineInsertPair(final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        if(!assembly1.junction().Chromosome.equals(assembly2.junction().Chromosome))
            return false;

        return withinLineProximity(
                assembly1.junction().Position, assembly2.junction().Position, assembly1.junction().Orient, assembly2.junction().Orient);
    }

    public static boolean withinLineProximity(final int pos1, final int pos2, final Orientation orient1, final Orientation orient2)
    {
        if(orient1 == orient2)
            return false;

        int posDiff = abs(pos2 - pos1);

        if(pos1 < pos2 == orient1.isForward())
            return posDiff <= LINE_INDEL_MAX_GAP;
        else
            return posDiff <= LINE_INDEL_MAX_OVERLAP;
    }

    public static int compareJunctions(
            final String chr1, final String chr2, final int pos1, final int pos2, final Orientation orient1, final Orientation orient2)
    {
        if(!chr1.equals(chr2))
        {
            int firstChrRank = HumanChromosome.chromosomeRank(chr1);
            int secondChrRank = HumanChromosome.chromosomeRank(chr2);

            return firstChrRank < secondChrRank ? -1 : 1;
        }

        if(pos1 == pos2)
        {
            if(orient1 == orient2)
                return 0;

            return orient1 .isForward() ? -1 : 1;
        }

        return pos1 < pos2 ? -1 : 1;
    }

    public static void checkStandardNucleotides(final SAMRecord read)
    {
        // simplifies usage downstream in assembly for other non-standard letters
        for(int i = 0; i < read.getReadBases().length; ++i)
        {
            if(!Nucleotides.isValidDnaBase(read.getReadBases()[i]))
            {
                read.getReadBases()[i] = DNA_N_BYTE;
                read.getBaseQualities()[i] = BaseQualAdjustment.BASE_QUAL_MINIMUM;
            }
        }
    }

    public static byte[] createByteArray(final int length, final byte value)
    {
        final byte[] array = new byte[length];

        for(int i = 0; i < array.length; ++i)
        {
            array[i] = value;
        }

        return array;
    }
}
