package com.hartwig.hmftools.common.variant.snpeff;

import static com.hartwig.hmftools.common.variant.snpeff.VariantAnnotationTest.createVariantAnnotationBuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Map;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.junit.Test;

public class CanonicalAnnotationSelectorTest {

    @Test
    public void testSelectCorrectAnnotation() {
        final Map<String, String> geneTranscriptMap = Maps.newHashMap();
        geneTranscriptMap.put("GENE", "TRANSCRIPT1");

        VariantAnnotation correct = createVariantAnnotationBuilder().gene("GENE").featureID("TRANSCRIPT1").build();
        VariantAnnotation incorrect = createVariantAnnotationBuilder().gene("GENE").featureID("TRANSCRIPT2").build();

        CanonicalAnnotationSelector victim = new CanonicalAnnotationSelector(geneTranscriptMap);
        Optional<VariantAnnotation> selected = victim.canonical("GENE", Lists.newArrayList(incorrect, correct));
        assertTrue(selected.isPresent());
        assertEquals(correct, selected.get());

        selected = victim.canonical("GENE", Lists.newArrayList(correct, incorrect));
        assertTrue(selected.isPresent());
        assertEquals(correct, selected.get());

        selected = victim.canonical("GENE2", Lists.newArrayList(correct, incorrect));
        assertFalse(selected.isPresent());
    }

}
