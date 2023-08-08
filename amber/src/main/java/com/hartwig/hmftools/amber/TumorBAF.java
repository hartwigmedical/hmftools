package com.hartwig.hmftools.amber;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

public class TumorBAF implements GenomePosition
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;

    public int NormalReadDepth;
    public int NormalRefSupport;
    public int NormalAltSupport;

    public int TumorReadDepth;
    public int TumorRefSupport;
    public int TumorAltSupport;
    public int TumorAltQuality;
    public int TumorIndelCount;

    public TumorBAF(final String chromosome, final int position, final String ref, final String alt)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;

        NormalReadDepth = 0;
        NormalRefSupport = 0;
        NormalAltSupport = 0;

        TumorReadDepth = 0;
        TumorRefSupport = 0;
        TumorAltSupport = 0;
        TumorAltQuality = 0;
        TumorIndelCount = 0;
    }

    public String chromosome() { return Chromosome; }
    public int position() { return Position; }

    public double refFrequency() { return TumorRefSupport / (double)TumorReadDepth; }
    public double altFrequency() {
        return TumorAltSupport / (double)TumorReadDepth;
    }
}
