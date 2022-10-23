package com.hartwig.hmftools.protect.algo;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectTestFactory;
import com.hartwig.serve.datamodel.ImmutableTreatment;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class EvidenceReportingCurationTest {

    @Test
    public void canApplyGeneBlacklist() {
        String gene1 = "any gene";
        String gene2 = "TP53";

        ProtectEvidence evidence1 = ProtectTestFactory.builder().gene(gene1).reported(true).build();
        ProtectEvidence evidence2 = ProtectTestFactory.builder().gene(gene2).reported(true).build();

        List<ProtectEvidence> evidence = EvidenceReportingCuration.applyReportingBlacklist(Lists.newArrayList(evidence1, evidence2));
        assertEquals(2, evidence.size());
        assertTrue(evidence.contains(evidence1));

        ProtectEvidence blacklisted = findByGene(evidence, gene2);
        assertFalse(blacklisted.reported());
    }

    @Test
    public void canApplyTreatmentBlacklist() {
        String treatment1 = "Chemotherapy";
        String treatment2 = "Immunotherapy";

        ProtectEvidence evidence1 = ProtectTestFactory.builder()
                .treatment(ImmutableTreatment.builder()
                        .treament(treatment1)
                        .sourceRelevantTreatmentApproaches(Sets.newHashSet("AA"))
                        .relevantTreatmentApproaches(Sets.newHashSet("A"))
                        .build())
                .reported(true)
                .build();
        ProtectEvidence evidence2 = ProtectTestFactory.builder()
                .treatment(ImmutableTreatment.builder()
                        .treament(treatment2)
                        .sourceRelevantTreatmentApproaches(Sets.newHashSet("AA"))
                        .relevantTreatmentApproaches(Sets.newHashSet("A"))
                        .build())
                .reported(true)
                .build();

        List<ProtectEvidence> evidence = EvidenceReportingCuration.applyReportingBlacklist(Lists.newArrayList(evidence1, evidence2));
        assertEquals(2, evidence.size());
        assertTrue(evidence.contains(evidence2));

        ProtectEvidence blacklisted = findByTreatment(evidence, treatment1);
        assertFalse(blacklisted.reported());
    }

    @NotNull
    private static ProtectEvidence findByTreatment(@NotNull Iterable<ProtectEvidence> evidences, @NotNull String treatment) {
        for (ProtectEvidence evidence : evidences) {
            if (evidence.treatment().treament().equals(treatment)) {
                return evidence;
            }
        }

        throw new IllegalStateException("Could not find evidence with treatment: " + treatment);
    }

    @NotNull
    private static ProtectEvidence findByGene(@NotNull Iterable<ProtectEvidence> evidences, @NotNull String gene) {
        for (ProtectEvidence evidence : evidences) {
            if (gene.equals(evidence.gene())) {
                return evidence;
            }
        }

        throw new IllegalStateException("Could not find evidence with gene: " + gene);
    }
}