package com.hartwig.hmftools.common.variant.impact;

import com.hartwig.hmftools.common.variant.CodingEffect;

public class VariantImpact
{
    public final String CanonicalGeneName;
    public final String CanonicalEffect;
    public final String CanonicalTranscript;
    public final CodingEffect CanonicalCodingEffect;
    public final String CanonicalHgvsCodingImpact;
    public final String CanonicalHgvsProteinImpact;
    public final boolean CanonicalSpliceRegion;

    public final String OtherReportableEffects;
    public final CodingEffect WorstCodingEffect;
    public final int GenesAffected;

    public VariantImpact(
            final String canonicalGeneName, final String canonicalTranscript, final String canonicalEffect,
            final CodingEffect canonicalCodingEffect, final String canonicalHgvsCodingImpact, final String canonicalHgvsProteinImpact,
            boolean canonicalSpliceRegion, final String otherReportableEffects, final CodingEffect worstCodingEffect, int genesAffected)
    {
        CanonicalGeneName = canonicalGeneName;
        CanonicalTranscript = canonicalTranscript;
        CanonicalEffect = canonicalEffect;
        CanonicalCodingEffect = canonicalCodingEffect;
        CanonicalSpliceRegion = canonicalSpliceRegion;
        CanonicalHgvsCodingImpact = canonicalHgvsCodingImpact;
        CanonicalHgvsProteinImpact = canonicalHgvsProteinImpact;
        OtherReportableEffects = otherReportableEffects;
        WorstCodingEffect = worstCodingEffect;
        GenesAffected = genesAffected;
    }

    public String gene()
    {
        return CanonicalGeneName;
    }

    public boolean equals(final VariantImpact other)
    {
        return CanonicalGeneName.equals(other.CanonicalGeneName) &&
                CanonicalTranscript.equals(other.CanonicalTranscript) &&
                CanonicalEffect.equals(other.CanonicalEffect) &&
                CanonicalCodingEffect == other.CanonicalCodingEffect &&
                CanonicalHgvsCodingImpact.equals(other.CanonicalHgvsCodingImpact) &&
                CanonicalHgvsProteinImpact.equals(other.CanonicalHgvsProteinImpact) &&
                CanonicalSpliceRegion == other.CanonicalSpliceRegion &&
                OtherReportableEffects.equals(other.OtherReportableEffects) &&
                WorstCodingEffect == other.WorstCodingEffect &&
                GenesAffected == other.GenesAffected;
    }

}
