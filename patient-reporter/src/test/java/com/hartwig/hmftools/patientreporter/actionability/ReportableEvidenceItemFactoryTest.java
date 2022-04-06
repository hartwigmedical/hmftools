package com.hartwig.hmftools.patientreporter.actionability;

import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ImmutableProtectSource;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceType;
import com.hartwig.hmftools.common.protect.ProtectTestFactory;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class ReportableEvidenceItemFactoryTest {

    @Test
    public void reportableFactoryWorksForTrivialCase() {
        ProtectEvidence item1 = ProtectTestFactory.testEvidenceBuilder()
                .event("A")
                .treatment("A")
                .onLabel(true)
                .level(EvidenceLevel.A)
                .protectSources(Sets.newHashSet(ImmutableProtectSource.builder()
                        .sources(Knowledgebase.VICC_CIVIC)
                        .sourceEvent(Strings.EMPTY)
                        .sourceUrls(Sets.newHashSet())
                        .evidenceType(ProtectEvidenceType.AMPLIFICATION)
                        .build()))
                .build();
        ProtectEvidence item2 = ProtectTestFactory.testEvidenceBuilder()
                .event("A")
                .treatment("A")
                .onLabel(true)
                .level(EvidenceLevel.A)
                .protectSources(Sets.newHashSet(ImmutableProtectSource.builder()
                        .sources(Knowledgebase.VICC_CGI)
                        .sourceEvent(Strings.EMPTY)
                        .sourceUrls(Sets.newHashSet())
                        .evidenceType(ProtectEvidenceType.AMPLIFICATION)
                        .build()))
                .build();
        ProtectEvidence item3 = ProtectTestFactory.testEvidenceBuilder()
                .event("B")
                .treatment("B")
                .onLabel(true)
                .level(EvidenceLevel.A)
                .protectSources(Sets.newHashSet(ImmutableProtectSource.builder()
                        .sources(Knowledgebase.ICLUSION)
                        .sourceEvent(Strings.EMPTY)
                        .sourceUrls(Sets.newHashSet())
                        .evidenceType(ProtectEvidenceType.AMPLIFICATION)
                        .build()))
                .build();
        ProtectEvidence item4 = ProtectTestFactory.testEvidenceBuilder()
                .event("C")
                .treatment("C")
                .onLabel(true)
                .level(EvidenceLevel.C)
                .protectSources(Sets.newHashSet(ImmutableProtectSource.builder()
                        .sources(Knowledgebase.VICC_CGI)
                        .sourceEvent(Strings.EMPTY)
                        .sourceUrls(Sets.newHashSet())
                        .evidenceType(ProtectEvidenceType.AMPLIFICATION)
                        .build()))
                .build();

        List<ProtectEvidence> nonTrials =
                ReportableEvidenceItemFactory.extractNonTrialsOnLabel(Lists.newArrayList(item1, item2, item3, item4));
        assertTrue(nonTrials.contains(item1));
        assertTrue(nonTrials.contains(item2));
        assertTrue(nonTrials.contains(item4));
    }
}