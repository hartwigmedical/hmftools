package com.hartwig.hmftools.pave.compare;

import static com.hartwig.hmftools.common.variant.VariantConsequence.EXON_LOSS_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.FUSION;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INITIATOR_CODON_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INTRAGENIC_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.NON_CANONICAL_START_CODON;
import static com.hartwig.hmftools.common.variant.VariantConsequence.PROTEIN_PROTEIN_CONTACT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.REGULATORY_REGION_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SEQUENCE_FEATURE;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SPLICE_REGION_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.STRUCTURAL_INTERACTION_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.TFBS_ABLATION;
import static com.hartwig.hmftools.common.variant.VariantConsequence.TRANSCRIPT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.TRANSCRIPT_ABLATION;
import static com.hartwig.hmftools.common.variant.VariantConsequence.TRANSCRIPT_AMPLIFICATION;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.pave.VariantTransImpact;

public final class ComparisonUtils
{
    public static List<SnpEffAnnotation> findMatchingAnnotations(final VariantTransImpact transImpact, final List<SnpEffAnnotation> annotations)
    {
        return annotations.stream()
                .filter(x -> x.featureID() != null && x.featureID().equals(transImpact.TransData.TransName))
                .filter(x -> !ignoreSnpEffAnnotation(x))
                .collect(Collectors.toList());
    }

    public static boolean ignoreSnpEffAnnotation(final SnpEffAnnotation annotation)
    {
        return annotation.consequences().size() == 1
                && (ignoreSnpEffConsequence(annotation.consequences().get(0)) || ignoreSnpEffEffect(annotation.effects()));
    }

    private static boolean ignoreSnpEffConsequence(final VariantConsequence consequence)
    {
        return consequence == SPLICE_REGION_VARIANT
            || consequence == SEQUENCE_FEATURE
            || consequence == TRANSCRIPT
            || consequence == REGULATORY_REGION_VARIANT
            || consequence == INITIATOR_CODON_VARIANT
            || consequence == EXON_LOSS_VARIANT
            || consequence == NON_CANONICAL_START_CODON
            || consequence == TRANSCRIPT_ABLATION
            || consequence == TFBS_ABLATION
            || consequence == STRUCTURAL_INTERACTION_VARIANT
            || consequence == FUSION
            || consequence == PROTEIN_PROTEIN_CONTACT
            || consequence == TRANSCRIPT_AMPLIFICATION;
    }

    private static boolean ignoreSnpEffEffect(final String effects)
    {
        return effects.equals("5_prime_UTR_premature_start_codon_gain_variant");
    }

    public static boolean effectsMatch(final VariantTransImpact transImpact, final SnpEffAnnotation annotation)
    {
        if(transImpact.effectsStr().equals(annotation.effects()))
            return true;

        // cannot compare every entry since some in SnpEff aren't valid in Pave
        List<VariantConsequence> consequences = annotation.consequences().stream()
                .filter(x -> !ignoreSnpEffConsequence(x))
                .collect(Collectors.toList());

        if(transImpact.effects().size() != consequences.size())
            return false;

        if(transImpact.effects().stream().noneMatch(x -> consequences.stream().anyMatch(y -> effectsEqual(x, y))))
            return false;

        return true;
    }

    public static boolean effectsEqual(final VariantEffect effect, final VariantConsequence consequence)
    {
        return consequence.isParentTypeOf(effect.effect());
    }

}
