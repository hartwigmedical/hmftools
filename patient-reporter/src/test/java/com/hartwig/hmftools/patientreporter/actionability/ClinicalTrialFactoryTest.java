package com.hartwig.hmftools.patientreporter.actionability;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ActionabilitySource;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.actionability.ImmutableEvidenceItem;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ClinicalTrialFactoryTest {

    @Test
    public void canCreateClinicalTrials() {
        EvidenceItem item = testEvidenceBuilder().event("event")
                .drug("acronym")
                .source(ActionabilitySource.ICLUSION)
                .reference("reference")
                .isOnLabel(false)
                .cancerType("Colorectal adenocarcinoma")
                .build();

        List<ClinicalTrial> trial = ClinicalTrialFactory.extractTrials(Lists.newArrayList(item));

        assertEquals(1, trial.size());
        assertEquals("event", trial.get(0).event());
        assertEquals("acronym", trial.get(0).acronym());
        assertEquals("reference", trial.get(0).reference());
        assertEquals(ActionabilitySource.ICLUSION, trial.get(0).source());
        assertFalse(trial.get(0).isOnLabel());
    }

    @NotNull
    private static ImmutableEvidenceItem.Builder testEvidenceBuilder() {
        return ImmutableEvidenceItem.builder()
                .source(ActionabilitySource.ICLUSION)
                .level(EvidenceLevel.LEVEL_A)
                .response(Strings.EMPTY)
                .drugsType(Strings.EMPTY);
    }
}