package com.hartwig.hmftools.isofox.data_loaders;

import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.util.Map;

public class SpliceVariant
{
    public final String GeneName;
    public final String Chromosome;

    public final int Position;
    public final String Ref;
    public final String Alt;
    public final String TransName;
    public final String CodingEffect;

    public SpliceVariant(
            final String geneName, final String chromosome, int position, final String ref, final String alt,
            final String transName, final String codingEffect)
    {
        GeneName = geneName;
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        TransName = transName;
        CodingEffect = codingEffect;
    }

    public static SpliceVariant fromCsv(final String[] items, final Map<String,Integer> fieldIndexMap)
    {
        return new SpliceVariant(
                items[fieldIndexMap.get("GeneName")],
                items[fieldIndexMap.get("Chromosome")],
                Integer.parseInt(items[fieldIndexMap.get("Position")]),
                items[fieldIndexMap.get("Ref")],
                items[fieldIndexMap.get("Alt")],
                items[fieldIndexMap.get("TransName")],
                items[fieldIndexMap.get("CodingEffect")]);

    }
}
