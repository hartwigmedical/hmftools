package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.CodingEffect.SPLICE;
import static com.hartwig.hmftools.common.variant.CodingEffect.SYNONYMOUS;
import static com.hartwig.hmftools.common.variant.CodingEffect.effect;
import static com.hartwig.hmftools.common.variant.VariantConsequence.FRAMESHIFT_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INFRAME_DELETION;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INFRAME_INSERTION;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INTRON_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.MISSENSE_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.PROTEIN_PROTEIN_CONTACT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SPLICE_ACCEPTOR_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SPLICE_DONOR_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SPLICE_REGION_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.STOP_GAINED;
import static com.hartwig.hmftools.common.variant.VariantConsequence.STRUCTURAL_INTERACTION_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SYNONYMOUS_VARIANT;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CodingEffectTest {

    private static final String RANDOM_GENE = "RANDOM";

    @Test
    public void testSingleEffect() {
        assertEffect(NONSENSE_OR_FRAMESHIFT, RANDOM_GENE, FRAMESHIFT_VARIANT);
        assertEffect(NONSENSE_OR_FRAMESHIFT, RANDOM_GENE, STOP_GAINED);

        assertEffect(MISSENSE, RANDOM_GENE, MISSENSE_VARIANT);
        assertEffect(MISSENSE, RANDOM_GENE, PROTEIN_PROTEIN_CONTACT);
        assertEffect(MISSENSE, RANDOM_GENE, STRUCTURAL_INTERACTION_VARIANT);
        assertEffect(MISSENSE, RANDOM_GENE, INFRAME_DELETION);
        assertEffect(MISSENSE, RANDOM_GENE, INFRAME_INSERTION);

        assertEffect(SPLICE, RANDOM_GENE, SPLICE_ACCEPTOR_VARIANT);
        assertEffect(SPLICE, RANDOM_GENE, SPLICE_DONOR_VARIANT);
        assertEffect(SYNONYMOUS, RANDOM_GENE, SYNONYMOUS_VARIANT);

        assertEffect(NONE, RANDOM_GENE, SPLICE_REGION_VARIANT);
        assertEffect(NONE, RANDOM_GENE, INTRON_VARIANT);
    }

    @Test
    public void testEffectPriority() {
        assertEffect(NONSENSE_OR_FRAMESHIFT, RANDOM_GENE, STOP_GAINED, MISSENSE_VARIANT, SPLICE_ACCEPTOR_VARIANT, INTRON_VARIANT);
        assertEffect(MISSENSE, RANDOM_GENE, MISSENSE_VARIANT, SPLICE_ACCEPTOR_VARIANT, SYNONYMOUS_VARIANT, INTRON_VARIANT);
        assertEffect(SPLICE, RANDOM_GENE, SPLICE_ACCEPTOR_VARIANT, INTRON_VARIANT);
        assertEffect(SYNONYMOUS, RANDOM_GENE, SYNONYMOUS_VARIANT, INTRON_VARIANT);
        assertEffect(NONE, RANDOM_GENE, INTRON_VARIANT);
    }

    @Test
    public void testTP53SpliceRegionVariant() {
        assertEffect(NONE, RANDOM_GENE, SPLICE_REGION_VARIANT);
        assertEffect(SPLICE, "TP53", SPLICE_REGION_VARIANT);
    }

    private static void assertEffect(@NotNull final CodingEffect expected, @NotNull final String gene, @NotNull final VariantConsequence... consequences) {
        assertEquals(expected, effect(gene, Lists.newArrayList(consequences)));
    }
}
