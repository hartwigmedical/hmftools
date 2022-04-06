package com.hartwig.hmftools.protect;

import static com.hartwig.hmftools.common.protect.ProtectTestFactory.testEvidenceBuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ImmutableProtectSource;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceType;
import com.hartwig.hmftools.common.protect.ProtectSource;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class EvidenceConsolidationTest {

    @Test
    public void canConsolidate() {
        String treatment1 = "treatment1";
        String url1 = "url1";
        String url3 = "url3";
        Knowledgebase source1 = Knowledgebase.VICC_CGI;
        Knowledgebase source2 = Knowledgebase.VICC_CIVIC;

        String treatment2 = "treatment2";

        ProtectEvidence evidence1 = testEvidenceBuilder().treatment(treatment1)
                .protectSources(Sets.newHashSet(ImmutableProtectSource.builder()
                        .sources(source1)
                        .sourceEvent("amp")
                        .addSourceUrls(url1)
                        .evidenceType(ProtectEvidenceType.AMPLIFICATION)
                        .build()))
                .build();
        ProtectEvidence evidence2 = testEvidenceBuilder().treatment(treatment1)
             .protectSources(Sets.newHashSet(ImmutableProtectSource.builder()
                        .sources(source2)
                        .sourceEvent("amp")
                        .addSourceUrls(url1)
                        .evidenceType(ProtectEvidenceType.DELETION)
                        .build()))
                .build();
        ProtectEvidence evidence3 = testEvidenceBuilder().treatment(treatment2)
                .protectSources(Sets.newHashSet(ImmutableProtectSource.builder()
                        .sources(source2)
                        .sourceEvent("amp")
                        .addSourceUrls(url3)
                        .evidenceType(ProtectEvidenceType.AMPLIFICATION)
                        .build()))
                .build();

        List<ProtectEvidence> consolidated = EvidenceConsolidation.consolidate(Lists.newArrayList(evidence1, evidence2, evidence3));

        assertEquals(2, consolidated.size());
        ProtectEvidence consolidatedEvidence1 = findByTreatment(consolidated, treatment1);
        assertEquals(2, consolidatedEvidence1.protectSources().size());

        ProtectSource protectSource1 = findBySource(consolidatedEvidence1.protectSources(), Knowledgebase.VICC_CGI);
        assertEquals(Knowledgebase.VICC_CGI, protectSource1.sources());
        assertEquals("amp", protectSource1.sourceEvent());
        assertEquals(Sets.newHashSet(url1), protectSource1.sourceUrls());
        assertEquals(ProtectEvidenceType.AMPLIFICATION, protectSource1.evidenceType());
        assertNull(protectSource1.rangeRank());

        ProtectSource protectSource2 = findBySource(consolidatedEvidence1.protectSources(), Knowledgebase.VICC_CIVIC);
        assertEquals(Knowledgebase.VICC_CIVIC, protectSource2.sources());
        assertEquals("amp", protectSource2.sourceEvent());
        assertEquals(Sets.newHashSet(url1), protectSource2.sourceUrls());
        assertEquals(ProtectEvidenceType.DELETION, protectSource2.evidenceType());
        assertNull(protectSource2.rangeRank());


        ProtectEvidence consolidatedEvidence2 = findByTreatment(consolidated, treatment2);
        assertEquals(1, consolidatedEvidence2.protectSources().size());

        ProtectSource protectSource3 = findBySource(consolidatedEvidence2.protectSources(), Knowledgebase.VICC_CIVIC);
        assertEquals(Knowledgebase.VICC_CIVIC, protectSource3.sources());
        assertEquals("amp", protectSource3.sourceEvent());
        assertEquals(Sets.newHashSet(url3), protectSource3.sourceUrls());
        assertEquals(ProtectEvidenceType.AMPLIFICATION, protectSource3.evidenceType());
        assertNull(protectSource3.rangeRank());
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

    @NotNull
    private static ProtectSource findBySource(@NotNull Set<ProtectSource> sources, @NotNull Knowledgebase source) {
        for (ProtectSource protectSource : sources) {
            if (protectSource.sources() == source) {
                return protectSource;
            }
        }

        throw new IllegalStateException("Could not find evidence with source: " + source);
    }
}