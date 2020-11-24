package com.hartwig.hmftools.isofox.novel.cohort;

import java.util.Map;

import com.hartwig.hmftools.common.variant.VariantType;

public class SpliceVariant
{
    public final String GeneName;
    public final String Chromosome;

    public final int Position;
    public final String Ref;
    public final String Alt;
    public final VariantType Type;
    public final String CodingEffect;
    public final String HgvsCodingImpact;
    public final String TriNucContext;

    public SpliceVariant(
            final String geneName, final String chromosome, int position, final VariantType type,
            final String ref, final String alt, final String codingEffect, final String hgvsCodingImpact, final String triNucContext)
    {
        GeneName = geneName;
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Type = type;
        HgvsCodingImpact = hgvsCodingImpact;
        TriNucContext = triNucContext;
        CodingEffect = codingEffect;
    }

    public static SpliceVariant fromCsv(final String[] items, final Map<String,Integer> fieldIndexMap)
    {
        return new SpliceVariant(
                items[fieldIndexMap.get("GeneName")],
                items[fieldIndexMap.get("Chromosome")],
                Integer.parseInt(items[fieldIndexMap.get("Position")]),
                VariantType.valueOf(items[fieldIndexMap.get("Type")]),
                items[fieldIndexMap.get("Ref")],
                items[fieldIndexMap.get("Alt")],
                items[fieldIndexMap.get("CodingEffect")],
                items[fieldIndexMap.get("HgvsCodingImpact")],
                items[fieldIndexMap.get("TriNucContext")]);

    }

    public String toString()
    {
        return String.format("gene(%s) loc(%s:%d) type(%s)", GeneName, Chromosome, Position, Type);
    }
}
