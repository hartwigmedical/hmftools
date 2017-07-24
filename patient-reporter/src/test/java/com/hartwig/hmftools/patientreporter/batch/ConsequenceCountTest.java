package com.hartwig.hmftools.patientreporter.batch;

import static com.hartwig.hmftools.common.variant.VariantAnnotationTest.createVariantAnnotationBuilder;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantConsequence;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ConsequenceCountTest {

    @Test
    public void canCountConsequences() {
        final List<SomaticVariant> variants = Lists.newArrayList();

        final SomaticVariant.Builder builder = new SomaticVariant.Builder();
        variants.add(builder.annotations(Lists.newArrayList(
                createVariantAnnotationBuilder(VariantConsequence.FRAMESHIFT_VARIANT).build(),
                createVariantAnnotationBuilder(VariantConsequence.INFRAME_DELETION).build())).build());
        variants.add(builder.annotations(Lists.newArrayList(
                createVariantAnnotationBuilder(VariantConsequence.INFRAME_DELETION).build())).build());

        final Map<VariantConsequence, Integer> counts = ConsequenceCount.count(variants);
        assertEquals(VariantConsequence.values().length, counts.size());
        assertEquals(0, counts.get(VariantConsequence.MISSENSE_VARIANT).intValue());
        assertEquals(1, counts.get(VariantConsequence.FRAMESHIFT_VARIANT).intValue());
        assertEquals(2, counts.get(VariantConsequence.INFRAME_DELETION).intValue());
    }

    @NotNull
    private static List<VariantConsequence> list(@NotNull final VariantConsequence consequence) {
        return Lists.newArrayList(consequence);
    }
}