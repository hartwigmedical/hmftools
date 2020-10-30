package com.hartwig.hmftools.linx.annotators;

import static java.lang.Math.abs;

public class IndelData
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;
    public final String Microhomology;
    public final int RepeatCount;
    public final double Ploidy;

    public final boolean IsDelete;
    public final boolean IsInsert;
    public final int Length;

    private static final int HIGH_REPEAT_COUNT = 4;

    public IndelData(
            final String chromosome, int position, final String ref, final String alt, final String microhomology,
            int repeatCount, double ploidy)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Microhomology = microhomology;
        RepeatCount = repeatCount;
        Ploidy = ploidy;

        int refLength = ref.length();
        int altLength = alt.length();

        if(alt.contains(","))
        {
            // take the first alt's length if there are more than 1
            final String[] altSplits = alt.split(",");
            altLength = altSplits[0].length();
        }

        IsDelete = altLength < refLength;
        IsInsert = !IsDelete;
        Length = abs(altLength - refLength);
    }

    public boolean highRepeatCount() { return RepeatCount >= HIGH_REPEAT_COUNT; }

    public boolean microhomology() { return IsDelete && highRepeatCount() && !Microhomology.isEmpty(); }

    public String toString()
    {
        return String.format("indel(%s %s:%d) repeat(%d) microhomology(%s) ploidy(%.2f)",
                IsDelete ? "DEL" : "INS", Chromosome, Position, RepeatCount, microhomology(), Ploidy);
    }

    public static final int CSV_REQUIRED_FIELDS = 8;
    public static final int INDEL_COL_SAMPLE = 0;

    public static IndelData fromString(final String[] indelData)
    {
        // parse CSV data
        int index = INDEL_COL_SAMPLE+1;

        return new IndelData(indelData[index++],
                Integer.parseInt(indelData[index++]),
                indelData[index++],
                indelData[index++],
                indelData[index++],
                Integer.parseInt(indelData[index++]),
                Double.parseDouble(indelData[index++]));
    }

}
