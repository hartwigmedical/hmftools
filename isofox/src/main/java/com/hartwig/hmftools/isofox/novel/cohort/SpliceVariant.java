package com.hartwig.hmftools.isofox.novel.cohort;

import java.util.Map;

import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.isofox.fusion.ChrGeneCollectionPair;

public class SpliceVariant
{
    public final String GeneName;
    public final String Chromosome;

    public final int Position;
    public final String Ref;
    public final String Alt;
    public final VariantType Type;
    public final String TransName;
    public final String CodingEffect;

    public SpliceVariant(
            final String geneName, final String chromosome, int position, final VariantType type,
            final String ref, final String alt, final String transName, final String codingEffect)
    {
        GeneName = geneName;
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Type = type;
        TransName = transName;
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
                items[fieldIndexMap.get("TransName")],
                items[fieldIndexMap.get("CodingEffect")]);

    }

    public String toString()
    {
        return String.format("gene(%s) loc(%s:%d) type(%s)", GeneName, Chromosome, Position, Type);
    }
}
