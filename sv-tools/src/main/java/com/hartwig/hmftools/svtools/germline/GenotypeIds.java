package com.hartwig.hmftools.svtools.germline;

public class GenotypeIds
{
    public final int ReferenceOrdinal;
    public final int TumorOrdinal;
    public final String ReferenceId;
    public final String TumorId;

    public GenotypeIds(final int referenceOrdinal, final int tumorOrdinal, final String referenceId, final String tumorId)
    {
        ReferenceOrdinal = referenceOrdinal;
        TumorOrdinal = tumorOrdinal;
        ReferenceId = referenceId;
        TumorId = tumorId;
    }

    public boolean hasReference() { return ReferenceOrdinal >= 0; }
}
