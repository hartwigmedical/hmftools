package com.hartwig.hmftools.common.variant.impact;

import static com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo.VAR_IMPACT_OTHER_REPORT_DELIM;
import static com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo.VAR_IMPACT_OTHER_REPORT_ITEM_COUNT;
import static com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo.VAR_IMPACT_OTHER_REPORT_ITEM_DELIM;

import java.util.List;

import com.hartwig.hmftools.common.variant.CodingEffect;

import org.apache.commons.compress.utils.Lists;

public class VariantImpact
{
    public final String CanonicalGeneName;
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
            final String canonicalGeneName, final String canonicalTranscript, final String canonicalEffect,
            final CodingEffect canonicalCodingEffect, final String canonicalHgvsCoding, final String canonicalHgvsProtein,
            boolean canonicalSpliceRegion, final String otherReportableEffects, final CodingEffect worstCodingEffect, int genesAffected)
    {
        CanonicalGeneName = canonicalGeneName;
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
        return CanonicalGeneName.equals(other.CanonicalGeneName) &&
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
                CanonicalGeneName, CanonicalCodingEffect, CanonicalEffect, WorstCodingEffect, CanonicalHgvsCoding, CanonicalHgvsProtein);
    }

}
