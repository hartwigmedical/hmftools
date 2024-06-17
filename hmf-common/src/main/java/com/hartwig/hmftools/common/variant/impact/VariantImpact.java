package com.hartwig.hmftools.common.variant.impact;

import com.hartwig.hmftools.common.variant.CodingEffect;

public class VariantImpact
{
    public final String GeneName;
    public final String CanonicalEffect;
    public final String CanonicalTranscript;
    public final CodingEffect CanonicalCodingEffect;
    public final String CanonicalHgvsCoding;
    public final String CanonicalHgvsProtein;
    public final boolean CanonicalSpliceRegion;

    public final String OtherReportableEffects;
    public final CodingEffect WorstCodingEffect;
    public final int GenesAffected;

    public VariantImpact(
            final String geneName, final String canonicalTranscript, final String canonicalEffect,
            final CodingEffect canonicalCodingEffect, final String canonicalHgvsCoding, final String canonicalHgvsProtein,
            boolean canonicalSpliceRegion, final String otherReportableEffects, final CodingEffect worstCodingEffect, int genesAffected)
    {
        GeneName = geneName;
        CanonicalTranscript = canonicalTranscript;
        CanonicalEffect = canonicalEffect;
        CanonicalCodingEffect = canonicalCodingEffect;
        CanonicalSpliceRegion = canonicalSpliceRegion;
        CanonicalHgvsCoding = canonicalHgvsCoding;
        CanonicalHgvsProtein = canonicalHgvsProtein;
        OtherReportableEffects = otherReportableEffects;
        WorstCodingEffect = worstCodingEffect;
        GenesAffected = genesAffected;
    }

    public boolean equals(final VariantImpact other)
    {
        return GeneName.equals(other.GeneName) &&
                CanonicalTranscript.equals(other.CanonicalTranscript) &&
                CanonicalEffect.equals(other.CanonicalEffect) &&
                CanonicalCodingEffect == other.CanonicalCodingEffect &&
                CanonicalHgvsCoding.equals(other.CanonicalHgvsCoding) &&
                CanonicalHgvsProtein.equals(other.CanonicalHgvsProtein) &&
                CanonicalSpliceRegion == other.CanonicalSpliceRegion &&
                OtherReportableEffects.equals(other.OtherReportableEffects) &&
                WorstCodingEffect == other.WorstCodingEffect &&
                GenesAffected == other.GenesAffected;
    }

    public String toString()
    {
        return String.format("%s canonical(%s: %s) worst(%s) hgvs(%s %s)",
                GeneName, CanonicalCodingEffect, CanonicalEffect, WorstCodingEffect, CanonicalHgvsCoding, CanonicalHgvsProtein);
    }

}
