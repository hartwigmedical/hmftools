package com.hartwig.hmftools.sage.vcf;

import static java.lang.String.format;

import java.util.StringJoiner;

import com.hartwig.hmftools.sage.common.VariantReadContext;

public class ReadContextVcfInfo
{
    public final int AlignmentStart;
    public final int VarIndex;
    public final String LeftFlank;
    public final String Core;
    public final String RightFlank;
    public final String Cigar;

    public static final String ITEM_DELIM = "-";

    public ReadContextVcfInfo(
            final int alignmentStart, final int varIndex, final String leftFlank, final String core,
            final String rightFlank, final String cigar)
    {
        AlignmentStart = alignmentStart;
        VarIndex = varIndex;
        LeftFlank = leftFlank;
        Core = core;
        RightFlank = rightFlank;
        Cigar = cigar;
    }

    public ReadContextVcfInfo(final VariantReadContext readContext)
    {
        AlignmentStart = readContext.AlignmentStart;
        VarIndex = readContext.VarIndex;
        LeftFlank = readContext.leftFlankStr();
        Core = readContext.coreStr();
        RightFlank = readContext.rightFlankStr();
        Cigar = readContext.readCigar();
    }

    public String readBases() { return LeftFlank + Core + RightFlank;}

    public String toVcfTag()
    {
        StringJoiner sj = new StringJoiner(ITEM_DELIM);
        sj.add(String.valueOf(AlignmentStart));
        sj.add(String.valueOf(VarIndex));
        sj.add(LeftFlank);
        sj.add(Core);
        sj.add(RightFlank);
        sj.add(Cigar);
        return sj.toString();
    }

    public String toString()
    {
        return format("index(%d) alignStart(%d) readBases(%s-%s-%s) cigar(%s)",
                VarIndex, AlignmentStart, LeftFlank, Core, RightFlank, Cigar);
    }

    public static ReadContextVcfInfo fromVcfTag(final String vcfTag)
    {
        String[] values = vcfTag.split(ITEM_DELIM, 6);

        if(values.length != 6)
            return null;

        return new ReadContextVcfInfo(
                Integer.parseInt(values[0]), Integer.parseInt(values[1]), values[2], values[3], values[4], values[5]);
    }
}
