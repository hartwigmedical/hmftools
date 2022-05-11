package com.hartwig.hmftools.protect;

import static com.hartwig.hmftools.common.protect.ProtectTestFactory.builder;

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
        Knowledgebase knowledgebase1 = Knowledgebase.VICC_CGI;
        Knowledgebase knowledgebase2 = Knowledgebase.VICC_CIVIC;

        String treatment2 = "treatment2";

        ProtectEvidence evidence1 = builder().treatment(treatment1)
                .sources(Sets.newHashSet(ImmutableProtectSource.builder()
                        .name(knowledgebase1)
                        .sourceEvent("amp")
                        .addSourceUrls(url1)
                        .evidenceType(ProtectEvidenceType.AMPLIFICATION)
                        .build()))
                .build();
        ProtectEvidence evidence2 = builder().treatment(treatment1)
             .sources(Sets.newHashSet(ImmutableProtectSource.builder()
                        .name(knowledgebase2)
                        .sourceEvent("amp")
                        .addSourceUrls(url1)
                        .evidenceType(ProtectEvidenceType.DELETION)
                        .build()))
                .build();
        ProtectEvidence evidence3 = builder().treatment(treatment2)
                .sources(Sets.newHashSet(ImmutableProtectSource.builder()
                        .name(knowledgebase2)
                        .sourceEvent("amp")
                        .addSourceUrls(url3)
                        .evidenceType(ProtectEvidenceType.AMPLIFICATION)
                        .build()))
                .build();

        List<ProtectEvidence> consolidated = EvidenceConsolidation.consolidate(Lists.newArrayList(evidence1, evidence2, evidence3));

        assertEquals(2, consolidated.size());
        ProtectEvidence consolidatedEvidence1 = findByTreatment(consolidated, treatment1);
        assertEquals(2, consolidatedEvidence1.sources().size());

        ProtectSource source1 = findByKnowledgebase(consolidatedEvidence1.sources(), Knowledgebase.VICC_CGI);
        assertEquals("amp", source1.sourceEvent());
        assertEquals(Sets.newHashSet(url1), source1.sourceUrls());
        assertEquals(ProtectEvidenceType.AMPLIFICATION, source1.evidenceType());
        assertNull(source1.rangeRank());

        ProtectSource source2 = findByKnowledgebase(consolidatedEvidence1.sources(), Knowledgebase.VICC_CIVIC);
        assertEquals("amp", source2.sourceEvent());
        assertEquals(Sets.newHashSet(url1), source2.sourceUrls());
        assertEquals(ProtectEvidenceType.DELETION, source2.evidenceType());
        assertNull(source2.rangeRank());

        ProtectEvidence consolidatedEvidence2 = findByTreatment(consolidated, treatment2);
        assertEquals(1, consolidatedEvidence2.sources().size());

        ProtectSource source3 = findByKnowledgebase(consolidatedEvidence2.sources(), Knowledgebase.VICC_CIVIC);
        assertEquals("amp", source3.sourceEvent());
        assertEquals(Sets.newHashSet(url3), source3.sourceUrls());
        assertEquals(ProtectEvidenceType.AMPLIFICATION, source3.evidenceType());
        assertNull(source3.rangeRank());
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
    private static ProtectSource findByKnowledgebase(@NotNull Set<ProtectSource> sources, @NotNull Knowledgebase knowledgebaseToFind) {
        for (ProtectSource source : sources) {
            if (source.name() == knowledgebaseToFind) {
                return source;
            }
        }

        throw new IllegalStateException("Could not find evidence from source: " + knowledgebaseToFind);
    }
}