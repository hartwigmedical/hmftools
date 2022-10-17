package com.hartwig.hmftools.gripss.common;

public class GenotypeIds
{
    public final int ReferenceOrdinal;
    public final int TumorOrdinal;
    public final String ReferenceId;
    public final String TumorId;
    public final boolean GermlineMode;

    public GenotypeIds(final int referenceOrdinal, final int tumorOrdinal, final String referenceId, final String tumorId, final boolean germlineMode)
    {
        ReferenceOrdinal = referenceOrdinal;
        TumorOrdinal = tumorOrdinal;
        ReferenceId = referenceId;
        TumorId = tumorId;
        GermlineMode = germlineMode;
    }

    public boolean hasReference() { return ReferenceOrdinal >= 0; }
}
