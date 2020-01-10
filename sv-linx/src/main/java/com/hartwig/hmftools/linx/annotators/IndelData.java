package com.hartwig.hmftools.linx.annotators;

import static java.lang.Math.abs;

public class IndelData
{
    public final String Chromosome;
    public final long Position;
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
            final String chromosome, long position, final String ref, final String alt, final String microhomology,
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

    public boolean microhology() { return IsDelete && highRepeatCount() && !Microhomology.isEmpty(); }

    public String toString()
    {
        return String.format("indel(%s %s:%d) repeat(%d) microhomology(%s) ploidy(%.2f)",
                IsDelete ? "DEL" : "INS", Chromosome, Position, RepeatCount, microhology(), Ploidy);
    }
}
