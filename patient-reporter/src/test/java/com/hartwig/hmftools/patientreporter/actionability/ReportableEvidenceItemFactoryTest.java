package com.hartwig.hmftools.patientreporter.actionability;

import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportableEvidenceItemFactoryTest {

    @Test
    public void reportableFactoryWorksForTrivialCase() {
        ProtectEvidence item1 = testProtectEvidenceBuilder().genomicEvent("A")
                .treatment("A")
                .level(com.hartwig.hmftools.common.serve.actionability.EvidenceLevel.A)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CIVIC))
                .build();
        ProtectEvidence item2 = testProtectEvidenceBuilder().genomicEvent("A")
                .treatment("A")
                .level(com.hartwig.hmftools.common.serve.actionability.EvidenceLevel.A)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CGI))
                .build();
        ProtectEvidence item3 = testProtectEvidenceBuilder().genomicEvent("B")
                .treatment("B")
                .level(com.hartwig.hmftools.common.serve.actionability.EvidenceLevel.A)
                .sources(Sets.newHashSet(Knowledgebase.ICLUSION))
                .build();
        ProtectEvidence item4 = testProtectEvidenceBuilder().genomicEvent("C")
                .treatment("C")
                .level(com.hartwig.hmftools.common.serve.actionability.EvidenceLevel.C)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CGI))
                .build();

        List<ProtectEvidence> nonTrials = ReportableEvidenceItemFactory.extractNonTrials(Lists.newArrayList(item1, item2, item3, item4));
        assertTrue(nonTrials.contains(item1));
        assertTrue(nonTrials.contains(item2));
        assertTrue(nonTrials.contains(item4));
    }

    @NotNull
    private static ImmutableProtectEvidence.Builder testProtectEvidenceBuilder() {
        return ImmutableProtectEvidence.builder()
                .germline(false)
                .reported(true)
                .onLabel(true)
                .direction(EvidenceDirection.RESPONSIVE)
                .urls(Sets.newHashSet("iclusion"));
    }
}