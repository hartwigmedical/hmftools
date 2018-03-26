package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.CodingEffect.SPLICE;
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

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CodingEffectTest {

    @Test
    public void testSingleEffect() {
        assertEffect(NONSENSE_OR_FRAMESHIFT, FRAMESHIFT_VARIANT);
        assertEffect(NONSENSE_OR_FRAMESHIFT, STOP_GAINED);

        assertEffect(MISSENSE, MISSENSE_VARIANT);
        assertEffect(MISSENSE, PROTEIN_PROTEIN_CONTACT);
        assertEffect(MISSENSE, STRUCTURAL_INTERACTION_VARIANT);
        assertEffect(MISSENSE, INFRAME_DELETION);
        assertEffect(MISSENSE, INFRAME_INSERTION);

        assertEffect(SPLICE, SPLICE_ACCEPTOR_VARIANT);
        assertEffect(SPLICE, SPLICE_DONOR_VARIANT);
        assertEffect(SPLICE, SPLICE_REGION_VARIANT);

        assertEffect(NONE, INTRON_VARIANT);
    }

    @Test
    public void testEffectPriority() {
        assertEffect(NONSENSE_OR_FRAMESHIFT, STOP_GAINED, SPLICE_ACCEPTOR_VARIANT, MISSENSE_VARIANT, INTRON_VARIANT);
        assertEffect(SPLICE, SPLICE_ACCEPTOR_VARIANT, MISSENSE_VARIANT, INTRON_VARIANT);
        assertEffect(MISSENSE, MISSENSE_VARIANT, INTRON_VARIANT);
        assertEffect(NONE, INTRON_VARIANT);
    }

    private void assertEffect(@NotNull CodingEffect expected, @NotNull final VariantConsequence... consequences) {
        assertEquals(expected, effect(Lists.newArrayList(consequences)));
    }
}
