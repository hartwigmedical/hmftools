package com.hartwig.hmftools.sage.impact;

import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.variant.VariantConsequence;

public class VariantImpact
{
    public final TranscriptData TransData;
    public final VariantConsequence Consequence;

    public VariantImpact(final TranscriptData transData, final VariantConsequence consequence)
    {
        TransData = transData;
        Consequence = consequence;
    }
}
