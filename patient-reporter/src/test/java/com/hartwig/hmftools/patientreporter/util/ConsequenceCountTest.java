package com.hartwig.hmftools.patientreporter.util;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantAnnotation;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.common.variant.VariantType;

import org.junit.Test;

public class ConsequenceCountTest {

    @Test
    public void canCountConsequences() {
        final List<SomaticVariant> variants = Lists.newArrayList();

        final SomaticVariant.Builder builder = new SomaticVariant.Builder(VariantType.SNP);
        variants.add(builder.annotations(Lists.newArrayList(
                new VariantAnnotation.Builder().consequence(VariantConsequence.FRAMESHIFT_VARIANT).build(),
                new VariantAnnotation.Builder().consequence(VariantConsequence.INFRAME_DELETION).build())).build());
        variants.add(builder.annotations(Lists.newArrayList(
                new VariantAnnotation.Builder().consequence(VariantConsequence.INFRAME_DELETION).build())).build());

        final Map<VariantConsequence, Integer> counts = ConsequenceCount.count(variants);
        assertEquals(VariantConsequence.values().length, counts.size());
        assertEquals(0, counts.get(VariantConsequence.INITIATOR_CODON_VARIANT).intValue());
        assertEquals(1, counts.get(VariantConsequence.FRAMESHIFT_VARIANT).intValue());
        assertEquals(2, counts.get(VariantConsequence.INFRAME_DELETION).intValue());
    }
}