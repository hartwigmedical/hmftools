package com.hartwig.hmftools.sage.impact;

import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public class ImpactClassifier
{
    private final RefGenomeInterface mRefGenome;

    public ImpactClassifier(final RefGenomeInterface refGenome)
    {
        mRefGenome = refGenome;
    }

    public VariantTransImpact classifyVariant(
            final TranscriptData transData)
    {
        return null;

    }

}
