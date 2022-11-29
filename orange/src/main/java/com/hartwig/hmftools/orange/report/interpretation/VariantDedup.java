package com.hartwig.hmftools.orange.report.interpretation;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.orange.algo.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.orange.algo.purple.PurpleVariant;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class VariantDedup {

    private static final Logger LOGGER = LogManager.getLogger(VariantDedup.class);

    private static final Set<VariantEffect> PHASED_EFFECTS =
            Sets.newHashSet(VariantEffect.PHASED_INFRAME_DELETION, VariantEffect.PHASED_INFRAME_INSERTION);

    private VariantDedup() {
    }

    @NotNull
    public static List<PurpleVariant> apply(@NotNull List<PurpleVariant> variants) {
        List<PurpleVariant> filtered = Lists.newArrayList();
        for (PurpleVariant variant : variants) {
            if (hasCanonicalPhasedEffect(variant) && hasSameEffectWithHigherVCN(variants, variant)) {
                LOGGER.debug("Dedup'ing variant '{}'", variant);
            } else {
                filtered.add(variant);
            }
        }
        return filtered;
    }

    private static boolean hasCanonicalPhasedEffect(@NotNull PurpleVariant variant) {
        for (VariantEffect effect : variant.canonicalImpact().effects()) {
            if (PHASED_EFFECTS.contains(effect)) {
                return true;
            }
        }
        return false;
    }

    private static boolean hasSameEffectWithHigherVCN(@NotNull List<PurpleVariant> variants, @NotNull PurpleVariant variantToMatch) {
        // We assume that variants with same effect have unique hgvs coding impact.
        Double minVariantCopyNumber = null;
        String uniqueHgvsCodingImpact = null;
        PurpleTranscriptImpact variantImpactToMatch = variantToMatch.canonicalImpact();

        for (PurpleVariant variant : variants) {
            PurpleTranscriptImpact variantImpact = variant.canonicalImpact();
            if (variantImpact.effects().equals(variantImpactToMatch.effects()) && variant.gene().equals(variantToMatch.gene())
                    && variantImpact.hgvsProteinImpact().equals(variantImpactToMatch.hgvsProteinImpact())) {
                if (minVariantCopyNumber == null || Doubles.lessThan(variant.variantCopyNumber(), minVariantCopyNumber)) {
                    minVariantCopyNumber = variant.variantCopyNumber();
                    uniqueHgvsCodingImpact = variantImpact.hgvsCodingImpact();
                } else if (Doubles.equal(variant.variantCopyNumber(), minVariantCopyNumber)) {
                    uniqueHgvsCodingImpact = variantImpact.hgvsCodingImpact().compareTo(uniqueHgvsCodingImpact) > 0
                            ? variantImpact.hgvsCodingImpact()
                            : uniqueHgvsCodingImpact;
                }
            }
        }

        boolean matchesMinAlleleCopyNumber = Doubles.equal(variantToMatch.variantCopyNumber(), minVariantCopyNumber);
        boolean matchesBestHgvsCodingImpact = variantImpactToMatch.hgvsCodingImpact().equals(uniqueHgvsCodingImpact);
        return !(matchesMinAlleleCopyNumber && matchesBestHgvsCodingImpact);
    }

}
