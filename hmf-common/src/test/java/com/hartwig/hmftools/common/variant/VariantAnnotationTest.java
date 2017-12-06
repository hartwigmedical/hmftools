package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class VariantAnnotationTest {

    @Test
    public void canGenerateConsequenceString() {
        final VariantAnnotation noConsequences = VariantAnnotationTest.createVariantAnnotationBuilder().build();
        assertEquals(Strings.EMPTY, noConsequences.consequenceString());

        final VariantAnnotation oneConsequence = createVariantAnnotationBuilder(VariantConsequence.MISSENSE_VARIANT).build();
        assertEquals(VariantConsequence.MISSENSE_VARIANT.readableSequenceOntologyTerm(), oneConsequence.consequenceString());

        final VariantAnnotation twoConsequences =
                createVariantAnnotationBuilder(VariantConsequence.MISSENSE_VARIANT, VariantConsequence.INFRAME_DELETION).build();
        assertTrue(twoConsequences.consequenceString().contains(VariantConsequence.MISSENSE_VARIANT.readableSequenceOntologyTerm()));
        assertTrue(twoConsequences.consequenceString().contains(VariantConsequence.INFRAME_DELETION.readableSequenceOntologyTerm()));

        final VariantAnnotation twoConsequencesIgnoreOne =
                createVariantAnnotationBuilder(VariantConsequence.MISSENSE_VARIANT, VariantConsequence.OTHER).build();
        assertEquals(VariantConsequence.MISSENSE_VARIANT.readableSequenceOntologyTerm(), twoConsequencesIgnoreOne.consequenceString());
    }

    @NotNull
    public static ImmutableVariantAnnotation.Builder createVariantAnnotationBuilder(@NotNull VariantConsequence... consequences) {
        return ImmutableVariantAnnotation.builder()
                .allele(Strings.EMPTY)
                .severity(Strings.EMPTY)
                .gene(Strings.EMPTY)
                .geneID(Strings.EMPTY)
                .featureType(Strings.EMPTY)
                .featureID(Strings.EMPTY)
                .transcriptBioType(Strings.EMPTY)
                .rank(Strings.EMPTY)
                .hgvsCoding(Strings.EMPTY)
                .hgvsProtein(Strings.EMPTY)
                .cDNAPosAndLength(Strings.EMPTY)
                .cdsPosAndLength(Strings.EMPTY)
                .aaPosAndLength(Strings.EMPTY)
                .distance(Strings.EMPTY)
                .addition(Strings.EMPTY)
                .consequences(Lists.newArrayList(consequences));
    }
}