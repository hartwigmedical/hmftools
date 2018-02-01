package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.SimplifiedEffect.INDEL_FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.SimplifiedEffect.INDEL_OTHER;
import static com.hartwig.hmftools.common.variant.SimplifiedEffect.NONE;
import static com.hartwig.hmftools.common.variant.SimplifiedEffect.effect;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class SimplifiedEffectTest {

    @Test
    public void testIndel() {
        assertEffect(INDEL_FRAMESHIFT, VariantType.INDEL, VariantConsequence.FRAMESHIFT_VARIANT);
        assertEffect(INDEL_FRAMESHIFT, VariantType.INDEL, VariantConsequence.FRAMESHIFT_VARIANT, VariantConsequence.PROTEIN_PROTEIN_CONTACT);
        assertEffect(INDEL_FRAMESHIFT, VariantType.INDEL, VariantConsequence.FRAMESHIFT_VARIANT, VariantConsequence.INTRON_VARIANT);
        assertEffect(INDEL_FRAMESHIFT, VariantType.INDEL, VariantConsequence.STOP_GAINED, VariantConsequence.INTRON_VARIANT);

        assertEffect(INDEL_OTHER, VariantType.INDEL, VariantConsequence.PROTEIN_PROTEIN_CONTACT);
        assertEffect(INDEL_OTHER, VariantType.INDEL, VariantConsequence.STRUCTURAL_INTERACTION_VARIANT, VariantConsequence.INTRON_VARIANT);
        assertEffect(INDEL_OTHER, VariantType.INDEL, VariantConsequence.INFRAME_DELETION, VariantConsequence.INTRON_VARIANT);
        assertEffect(INDEL_OTHER, VariantType.INDEL, VariantConsequence.INFRAME_INSERTION, VariantConsequence.INTRON_VARIANT);

        assertEffect(NONE, VariantType.INDEL, VariantConsequence.INTRON_VARIANT);
    }


    private void assertEffect(@NotNull SimplifiedEffect expected, @NotNull VariantType type, @NotNull final VariantConsequence... consequences) {
        assertEquals(expected, effect(type, Lists.newArrayList(consequences)));
    }

    @NotNull
    private static List<VariantConsequence> create(@NotNull final VariantConsequence... consequences) {
        return Lists.newArrayList(consequences);
    }
}
