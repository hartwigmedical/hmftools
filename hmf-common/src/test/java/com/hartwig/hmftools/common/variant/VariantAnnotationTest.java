package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class VariantAnnotationTest {

    @Test
    public void canGenerateConsequenceString() {
        final VariantAnnotation noConsequences = new VariantAnnotation.Builder().build();
        assertEquals(Strings.EMPTY, noConsequences.consequenceString());

        final VariantAnnotation oneConsequence = new VariantAnnotation.Builder().consequences(
                Lists.newArrayList(VariantConsequence.MISSENSE_VARIANT)).build();
        assertEquals(VariantConsequence.MISSENSE_VARIANT.readableSequenceOntologyTerm(),
                oneConsequence.consequenceString());

        final VariantAnnotation twoConsequences = new VariantAnnotation.Builder().consequences(
                Lists.newArrayList(VariantConsequence.MISSENSE_VARIANT, VariantConsequence.INFRAME_DELETION)).build();
        assertTrue(twoConsequences.consequenceString().contains(
                VariantConsequence.MISSENSE_VARIANT.readableSequenceOntologyTerm()));
        assertTrue(twoConsequences.consequenceString().contains(
                VariantConsequence.INFRAME_DELETION.readableSequenceOntologyTerm()));

        final VariantAnnotation twoConsequencesIgnoreOne = new VariantAnnotation.Builder().consequences(
                Lists.newArrayList(VariantConsequence.MISSENSE_VARIANT, VariantConsequence.OTHER)).build();
        assertEquals(VariantConsequence.MISSENSE_VARIANT.readableSequenceOntologyTerm(),
                twoConsequencesIgnoreOne.consequenceString());
    }
}