package com.hartwig.hmftools.pave;

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
}
