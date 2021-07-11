package com.hartwig.hmftools.sage.impact;

import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.variant.VariantConsequence;

public class VariantTransImpact
{
    public final TranscriptData TransData;
    public final VariantConsequence Consequence;

    public VariantTransImpact(final TranscriptData transData, final VariantConsequence consequence)
    {
        TransData = transData;
        Consequence = consequence;
    }

    public String codingChange() { return ""; }
    public String proteinChange() { return ""; }
}
