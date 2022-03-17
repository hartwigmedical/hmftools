package com.hartwig.hmftools.protect;

import static com.hartwig.hmftools.common.protect.ProtectTestFactory.testEvidenceBuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class EvidenceConsolidationTest {

    @Test
    public void canConsolidate() {
        String treatment1 = "treatment1";
        String url1 = "url1";
        String url2 = "url2";
        String url3 = "url3";
        Knowledgebase source1 = Knowledgebase.VICC_CGI;
        Knowledgebase source2 = Knowledgebase.VICC_CIVIC;

        String treatment2 = "treatment2";

        ProtectEvidence evidence1 = testEvidenceBuilder().treatment(treatment1).addSourceUrls(url1).addSources(source1).build();
        ProtectEvidence evidence2 = testEvidenceBuilder().treatment(treatment1).addSourceUrls(url2).addSources(source2).build();
        ProtectEvidence evidence3 = testEvidenceBuilder().treatment(treatment2).addSourceUrls(url3).addSources(source1).build();

        List<ProtectEvidence> consolidated = EvidenceConsolidation.consolidate(Lists.newArrayList(evidence1, evidence2, evidence3));

        assertEquals(2, consolidated.size());
        ProtectEvidence consolidatedEvidence1 = findByTreatment(consolidated, treatment1);
        assertTrue(consolidatedEvidence1.sourceUrls().contains(url1));
        assertTrue(consolidatedEvidence1.sourceUrls().contains(url2));
        assertFalse(consolidatedEvidence1.sourceUrls().contains(url3));
        assertTrue(consolidatedEvidence1.sources().contains(source1));
        assertTrue(consolidatedEvidence1.sources().contains(source2));

        ProtectEvidence consolidatedEvidence2 = findByTreatment(consolidated, treatment2);
        assertFalse(consolidatedEvidence2.sourceUrls().contains(url1));
        assertFalse(consolidatedEvidence2.sourceUrls().contains(url2));
        assertTrue(consolidatedEvidence2.sourceUrls().contains(url3));
        assertTrue(consolidatedEvidence2.sources().contains(source1));
        assertFalse(consolidatedEvidence2.sources().contains(source2));
    }

    @NotNull
    private static ProtectEvidence findByTreatment(@NotNull Iterable<ProtectEvidence> evidences, @NotNull String treatment) {
        for (ProtectEvidence evidence : evidences) {
            if (evidence.treatment().equals(treatment)) {
                return evidence;
            }
        }

        throw new IllegalStateException("Could not find evidence with treatment: " + treatment);
    }
}