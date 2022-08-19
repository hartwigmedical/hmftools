package com.hartwig.hmftools.pave.compare;

import static com.hartwig.hmftools.common.variant.CodingEffect.NONE;
import static com.hartwig.hmftools.common.variant.CodingEffect.UNDEFINED;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.SYNONYMOUS;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.UPSTREAM_GENE;
import static com.hartwig.hmftools.pave.HgvsCoding.HGVS_TYPE_DUP;
import static com.hartwig.hmftools.pave.HgvsCoding.HGVS_TYPE_INS;


import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.pave.VariantData;
import com.hartwig.hmftools.pave.VariantTransImpact;

public final class ComparisonUtils
{
   public static boolean hasCodingEffect(final CodingEffect effect)
    {
        return effect != UNDEFINED && effect != NONE;
    }

    public static boolean hasCodingEffectDiff(final VariantImpact variantImpact, final RefVariantData refVariant)
    {
        if(!hasCodingEffect(variantImpact.CanonicalCodingEffect) && !hasCodingEffect(refVariant.CanonicalCodingEffect))
            return false;

        if(refVariant.CanonicalCodingEffect == UNDEFINED)
            return false;

        return variantImpact.CanonicalCodingEffect != refVariant.CanonicalCodingEffect;
    }

    public static boolean hasHgvsCodingDiff(
            final VariantTransImpact transImpact, final VariantTransImpact raTransImpact, final RefVariantData refVariant)
    {
        if(refVariant.CanonicalCodingEffect == UNDEFINED) // different gene definitions
            return false;

        if(transImpact.hasEffect(UPSTREAM_GENE))
            return false;

        if(transImpact.hgvsCoding().equals(refVariant.HgvsCodingImpact))
            return false;

        if(raTransImpact != null && raTransImpact.hgvsCoding().equals(refVariant.HgvsCodingImpact))
            return false;

        if(transImpact.hgvsCoding().contains(HGVS_TYPE_DUP) && refVariant.HgvsCodingImpact.contains(HGVS_TYPE_INS))
            return false;

        return true;
    }

    public static boolean hasHgvsProteinDiff(
            final VariantData variant, final VariantTransImpact transImpact, final VariantTransImpact raTransImpact, final RefVariantData refVariant)
    {
        if(refVariant.CanonicalCodingEffect == UNDEFINED) // different gene definitions
            return false;

        // ignore synonymous since SnpEff doesn't use the 'equals' sign
        if(transImpact.hasEffect(SYNONYMOUS) && refVariant.CanonicalEffect.contains(SYNONYMOUS.effect()))
            return false;

        if(transImpact.hgvsProtein().equals(refVariant.HgvsProteinImpact))
            return false;

        if(raTransImpact != null && raTransImpact.hgvsProtein().equals(refVariant.HgvsProteinImpact))
            return false;

        // SnpEff uses a different convention for MNVs spanning 2 codons
        if(variant.isBaseChange() && transImpact.hasProteinContext() && transImpact.proteinContext().RefAminoAcids.length() > 1)
            return false;

        return true;
    }

}
