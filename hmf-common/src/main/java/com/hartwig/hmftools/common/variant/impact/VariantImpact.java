package com.hartwig.hmftools.common.variant.impact;

import com.hartwig.hmftools.common.variant.CodingEffect;

public class VariantImpact
{
    public final int GenesAffected;

    public final String CanonicalGeneName;
    public final String CanonicalGeneId;
    public final String CanonicalEffect;
    public final String CanonicalTranscript;
    public final CodingEffect CanonicalCodingEffect;
    public final String CanonicalHgvsCodingImpact;
    public final String CanonicalHgvsProteinImpact;

    public final String WorstGene;
    public final String WorstEffect;
    public final String WorstTranscript;
    public final CodingEffect WorstCodingEffect;

    public VariantImpact(
            final int genesAffected, final String canonicalGeneId, final String canonicalGeneName, final String canonicalEffect,
            final String canonicalTranscript, final CodingEffect canonicalCodingEffect, final String canonicalHgvsCodingImpact,
            final String canonicalHgvsProteinImpact, final String worstGene,
            final String worstEffect, final String worstTranscript, final CodingEffect worstCodingEffect)
    {
        GenesAffected = genesAffected;
        CanonicalGeneName = canonicalGeneName;
        CanonicalGeneId = canonicalGeneId;
        CanonicalEffect = canonicalEffect;
        CanonicalTranscript = canonicalTranscript;
        CanonicalCodingEffect = canonicalCodingEffect;
        CanonicalHgvsCodingImpact = canonicalHgvsCodingImpact;
        CanonicalHgvsProteinImpact = canonicalHgvsProteinImpact;
        WorstGene = worstGene;
        WorstEffect = worstEffect;
        WorstTranscript = worstTranscript;
        WorstCodingEffect = worstCodingEffect;
    }

    public String gene()
    {
        return CanonicalGeneName.isEmpty() ? WorstGene : CanonicalGeneName;
    }

}
