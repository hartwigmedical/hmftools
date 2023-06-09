package com.hartwig.hmftools.pave.annotation;

class PonVariantData
{
    public final String Ref;
    public final String Alt;
    public final int Samples;
    public final int MaxSampleReads;
    public final int TotalSampleReads;

    public PonVariantData(final String ref, final String alt, final int samples, final int maxSampleReads, final int totalSampleReads)
    {
        Ref = ref;
        Alt = alt;
        Samples = samples;
        MaxSampleReads = maxSampleReads;
        TotalSampleReads = totalSampleReads;
    }

    public boolean matches(final String ref, final String alt) { return ref.equals(Ref) && alt.equals(Alt); }

    public String toString()
    {
        return String.format("%s>%s sample(%d) maxReads(%d) total(%d)", Ref, Alt, Samples, MaxSampleReads, TotalSampleReads);
    }
}
