package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_INDEL_MAX_GAP;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_INDEL_MAX_OVERLAP;
import static com.hartwig.hmftools.esvee.common.SvConstants.LOW_BASE_QUAL_THRESHOLD;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

public final class CommonUtils
{
    public static boolean aboveMinQual(byte qual) { return qual >= LOW_BASE_QUAL_THRESHOLD; }
    public static boolean aboveMinQual(int qual) { return qual >= LOW_BASE_QUAL_THRESHOLD; }

    public static boolean belowMinQual(byte qual) { return qual < LOW_BASE_QUAL_THRESHOLD; }
    public static boolean belowMinQual(int qual) { return qual < LOW_BASE_QUAL_THRESHOLD; }

    public static boolean isDiscordantFragment(
            final SAMRecord read, final int fragmentLengthUpperBound, @Nullable final SupplementaryReadData suppData)
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

        int fragmentSize = abs(read.getInferredInsertSize());

        return fragmentSize == 0 || (fragmentLengthUpperBound > 0 && fragmentSize >= fragmentLengthUpperBound);
    }

    public static boolean isLineInsertPair(final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        if(!assembly1.junction().Chromosome.equals(assembly2.junction().Chromosome))
            return false;

        return withLineProximity(
                assembly1.junction().Position, assembly2.junction().Position, assembly1.junction().Orient, assembly2.junction().Orient);
    }

    public static boolean withLineProximity(final int pos1, final int pos2, final Orientation orient1, final Orientation orient2)
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

    public static StructuralVariantType formSvType(
            final String chrStart, final String chrEnd, final int posStart, final int posEnd,
            final Orientation orientStart, final Orientation orientEnd, final boolean hasInsertedBases)
    {
        if(!chrStart.equals(chrEnd))
            return BND;

        if(orientStart != orientEnd)
        {
            int posDiff = abs(posStart - posEnd);

            if(posDiff == 1 && hasInsertedBases)
                return INS;

            if(posDiff == 0)
                return DUP;

            boolean firstIsLower = posStart < posEnd;

            return (firstIsLower == orientStart.isForward()) ? DEL : DUP;
        }
        else
        {
            return INV;
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
