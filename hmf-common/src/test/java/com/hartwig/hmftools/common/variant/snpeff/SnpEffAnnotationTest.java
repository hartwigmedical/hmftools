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
        SnpEffAnnotation noConsequences = withConsequences();
        assertEquals(Strings.EMPTY, noConsequences.consequenceString());

        SnpEffAnnotation oneConsequence = withConsequences(VariantConsequence.MISSENSE_VARIANT);
        assertEquals(VariantConsequence.MISSENSE_VARIANT.readableSequenceOntologyTerm(), oneConsequence.consequenceString());

        SnpEffAnnotation twoConsequences = withConsequences(VariantConsequence.MISSENSE_VARIANT, VariantConsequence.INFRAME_DELETION);
        assertTrue(twoConsequences.consequenceString().contains(VariantConsequence.MISSENSE_VARIANT.readableSequenceOntologyTerm()));
        assertTrue(twoConsequences.consequenceString().contains(VariantConsequence.INFRAME_DELETION.readableSequenceOntologyTerm()));

        SnpEffAnnotation twoConsequencesIgnoreOne = withConsequences(VariantConsequence.MISSENSE_VARIANT, VariantConsequence.OTHER);
        assertEquals(VariantConsequence.MISSENSE_VARIANT.readableSequenceOntologyTerm(), twoConsequencesIgnoreOne.consequenceString());
    }

    @Test
    public void canStripVersionFromTranscript() {
        SnpEffAnnotation normal = withTranscript("ENST001");
        assertEquals("ENST001", normal.transcript());

        SnpEffAnnotation versioned = withTranscript("ENST002.12");
        assertEquals("ENST002", versioned.transcript());
    }

    @NotNull
    private static SnpEffAnnotation withConsequences(@NotNull VariantConsequence... consequences) {
        return SnpEffAnnotationTestFactory.builder().consequences(Lists.newArrayList(consequences)).build();
    }

    @NotNull
    private static SnpEffAnnotation withTranscript(@NotNull String transcript) {
        return SnpEffAnnotationTestFactory.builder().featureID(transcript).build();
    }
}