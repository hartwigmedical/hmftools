package com.hartwig.hmftools.common.variant.snpeff;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.VariantConsequence;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class SnpEffAnnotationTest {

    @Test
    public void canGenerateConsequenceString() {
        final SnpEffAnnotation noConsequences = annotation();
        assertEquals(Strings.EMPTY, noConsequences.consequenceString());

        final SnpEffAnnotation oneConsequence = annotation(VariantConsequence.MISSENSE_VARIANT);
        assertEquals(VariantConsequence.MISSENSE_VARIANT.readableSequenceOntologyTerm(), oneConsequence.consequenceString());

        final SnpEffAnnotation twoConsequences = annotation(VariantConsequence.MISSENSE_VARIANT, VariantConsequence.INFRAME_DELETION);
        assertTrue(twoConsequences.consequenceString().contains(VariantConsequence.MISSENSE_VARIANT.readableSequenceOntologyTerm()));
        assertTrue(twoConsequences.consequenceString().contains(VariantConsequence.INFRAME_DELETION.readableSequenceOntologyTerm()));

        final SnpEffAnnotation twoConsequencesIgnoreOne = annotation(VariantConsequence.MISSENSE_VARIANT, VariantConsequence.OTHER);
        assertEquals(VariantConsequence.MISSENSE_VARIANT.readableSequenceOntologyTerm(), twoConsequencesIgnoreOne.consequenceString());
    }

    @NotNull
    private static SnpEffAnnotation annotation(@NotNull VariantConsequence... consequences) {
        return SnpEffAnnotationTestFactory.builder().consequences(Lists.newArrayList(consequences)).build();
    }
}