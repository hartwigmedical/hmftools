package com.hartwig.hmftools.patientreporter.actionability;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ImmutableKnowledgebaseSource;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.EvidenceType;
import com.hartwig.hmftools.common.protect.ProtectTestFactory;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.serve.actionability.ImmutableTreatment;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class ClinicalTrialFactoryTest {

    @Test
    public void canExtractClinicalTrials() {
        ProtectEvidence evidence = ProtectTestFactory.builder()
                .event("event")
                .germline(false)
                .reported(true)
                .treatment(ImmutableTreatment.builder()
                        .treament("acronym")
                        .sourceRelevantTreatmentApproaches(Sets.newHashSet())
                        .relevantTreatmentApproaches(Sets.newHashSet())
                        .build())
                .onLabel(true)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(ImmutableKnowledgebaseSource.builder()
                        .name(Knowledgebase.ICLUSION)
                        .sourceEvent(Strings.EMPTY)
                        .sourceUrls(Sets.newHashSet())
                        .evidenceType(EvidenceType.AMPLIFICATION)
                        .build()))
                .build();

        List<ProtectEvidence> trial = ClinicalTrialFactory.extractOnLabelTrials(Lists.newArrayList(evidence));

        assertEquals(1, trial.size());
        assertEquals("event", trial.get(0).event());
        assertEquals("acronym", trial.get(0).treatment().treament());
    }
}