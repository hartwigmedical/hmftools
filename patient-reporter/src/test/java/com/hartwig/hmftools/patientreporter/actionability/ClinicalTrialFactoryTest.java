package com.hartwig.hmftools.patientreporter.actionability;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectTestFactory;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.junit.Test;

public class ClinicalTrialFactoryTest {

    @Test
    public void canExtractClinicalTrials() {
        ProtectEvidence evidence = ProtectTestFactory.testEvidenceBuilder()
                .event("event")
                .germline(false)
                .reported(true)
                .treatment("acronym")
                .onLabel(true)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.ICLUSION))
                .build();

        List<ProtectEvidence> trial = ClinicalTrialFactory.extractOnLabelTrials(Lists.newArrayList(evidence));

        assertEquals(1, trial.size());
        assertEquals("event", trial.get(0).genomicEvent());
        assertEquals("acronym", trial.get(0).treatment());
    }
}