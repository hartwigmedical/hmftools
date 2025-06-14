package com.hartwig.hmftools.common.sv;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;

import com.hartwig.hmftools.common.genome.region.Orientation;

import htsjdk.samtools.SAMRecord;

public final class SvUtils
{
    public static final int DEFAULT_DISCORDANT_FRAGMENT_LENGTH = 1000;

    // must match the small DEl-DUP threshold in Esvee
    public static final int SHORT_INDEL_LENGTH = 1000;

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

    public static boolean isIndel(final StructuralVariantType type)
    {
        return type == DEL || type == DUP || type == INS;
    }

    public static boolean hasShortIndelLength(final int svLength) { return svLength <= SHORT_INDEL_LENGTH; }

    public static boolean isShortLocalDelDupIns(final StructuralVariantType svType, final int svLength)
    {
        return isIndel(svType) && hasShortIndelLength(svLength);
    }

    public static boolean isDiscordant(final SAMRecord record) { return isDiscordant(record, DEFAULT_DISCORDANT_FRAGMENT_LENGTH); }

    public static boolean isDiscordant(final SAMRecord record, final int discordantPairFragmentLength)
    {
        if(!record.getReadPairedFlag())
            return false;

        if(!record.getReferenceName().equals(record.getMateReferenceName()))
            return true;

        if(record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag())
            return true;

        int fragmentSize = abs(record.getInferredInsertSize());

        return fragmentSize == 0 || fragmentSize >= discordantPairFragmentLength;
    }

}
