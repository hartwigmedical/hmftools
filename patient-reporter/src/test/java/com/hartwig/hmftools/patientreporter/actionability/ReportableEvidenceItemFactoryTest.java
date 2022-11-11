package com.hartwig.hmftools.patientreporter.actionability;

import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.EvidenceType;
import com.hartwig.hmftools.common.protect.ImmutableKnowledgebaseSource;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectTestFactory;
import com.hartwig.serve.datamodel.EvidenceLevel;
import com.hartwig.serve.datamodel.ImmutableTreatment;
import com.hartwig.serve.datamodel.Knowledgebase;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class ReportableEvidenceItemFactoryTest {

    @Test
    public void reportableFactoryWorksForTrivialCase() {
        ProtectEvidence item1 = ProtectTestFactory.builder()
                .event("A")
                .treatment(ImmutableTreatment.builder()
                        .treament("A")
                        .sourceRelevantTreatmentApproaches(Sets.newHashSet())
                        .relevantTreatmentApproaches(Sets.newHashSet())
                        .build())
                .onLabel(true)
                .level(EvidenceLevel.A)
                .sources(Sets.newHashSet(ImmutableKnowledgebaseSource.builder()
                        .name(Knowledgebase.VICC_CIVIC)
                        .sourceEvent(Strings.EMPTY)
                        .sourceUrls(Sets.newHashSet())
                        .evidenceType(EvidenceType.AMPLIFICATION)
                        .build()))
                .build();
        ProtectEvidence item2 = ProtectTestFactory.builder()
                .event("A")
                .treatment(ImmutableTreatment.builder()
                        .treament("A")
                        .sourceRelevantTreatmentApproaches(Sets.newHashSet())
                        .relevantTreatmentApproaches(Sets.newHashSet())
                        .build())
                .onLabel(true)
                .level(EvidenceLevel.A)
                .sources(Sets.newHashSet(ImmutableKnowledgebaseSource.builder()
                        .name(Knowledgebase.VICC_CGI)
                        .sourceEvent(Strings.EMPTY)
                        .sourceUrls(Sets.newHashSet())
                        .evidenceType(EvidenceType.AMPLIFICATION)
                        .build()))
                .build();
        ProtectEvidence item3 = ProtectTestFactory.builder()
                .event("B")
                .treatment(ImmutableTreatment.builder()
                        .treament("B")
                        .sourceRelevantTreatmentApproaches(Sets.newHashSet())
                        .relevantTreatmentApproaches(Sets.newHashSet())
                        .build())
                .onLabel(true)
                .level(EvidenceLevel.A)
                .sources(Sets.newHashSet(ImmutableKnowledgebaseSource.builder()
                        .name(Knowledgebase.ICLUSION)
                        .sourceEvent(Strings.EMPTY)
                        .sourceUrls(Sets.newHashSet())
                        .evidenceType(EvidenceType.AMPLIFICATION)
                        .build()))
                .build();
        ProtectEvidence item4 = ProtectTestFactory.builder()
                .event("C")
                .treatment(ImmutableTreatment.builder()
                        .treament("C")
                        .sourceRelevantTreatmentApproaches(Sets.newHashSet())
                        .relevantTreatmentApproaches(Sets.newHashSet())
                        .build())
                .onLabel(true)
                .level(EvidenceLevel.C)
                .sources(Sets.newHashSet(ImmutableKnowledgebaseSource.builder()
                        .name(Knowledgebase.VICC_CGI)
                        .sourceEvent(Strings.EMPTY)
                        .sourceUrls(Sets.newHashSet())
                        .evidenceType(EvidenceType.AMPLIFICATION)
                        .build()))
                .build();

        List<ProtectEvidence> nonTrials =
                ReportableEvidenceItemFactory.extractNonTrialsOnLabel(Lists.newArrayList(item1, item2, item3, item4));
        assertTrue(nonTrials.contains(item1));
        assertTrue(nonTrials.contains(item2));
        assertTrue(nonTrials.contains(item4));
    }
}