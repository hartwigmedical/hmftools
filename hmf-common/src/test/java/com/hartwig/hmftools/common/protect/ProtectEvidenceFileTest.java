package com.hartwig.hmftools.common.protect;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ProtectEvidenceFileTest {

    private static final String EVIDENCE_TSV = Resources.getResource("protect/example.tsv").getPath();

    @Test
    public void canReadProtectEvidenceFile() throws IOException {
        List<ProtectEvidence> evidences = ProtectEvidenceFile.read(EVIDENCE_TSV);
        assertEquals(5, evidences.size());

        // evidences 1
        ProtectEvidence evidence1 = findByTreatmentAndEvent(evidences, "Cobimetinib + Vemurafenib", "p.Val600Glu");
        assertEquals(1, evidence1.sources().size());

        ProtectSource evidence1Source = findBySource(evidence1.sources(), Knowledgebase.VICC_CGI);
        assertEquals("hotspot", evidence1Source.sourceEvent());
        assertEquals(Sets.newHashSet("https://www.google.com/#q=FDA"), evidence1Source.sourceUrls());
        assertEquals(ProtectEvidenceType.HOTSPOT_MUTATION, evidence1Source.evidenceType());
        assertNull(evidence1Source.rangeRank());
        assertTrue(evidence1Source.evidenceUrls().isEmpty());

        // evidences 2
        ProtectEvidence evidence2 = findByTreatmentAndEvent(evidences, "Dabrafenib", "p.Val600Glu");
        assertEquals(1, evidence2.sources().size());

        ProtectSource evidence2Source = findBySource(evidence2.sources(), Knowledgebase.VICC_CGI);
        assertEquals("hotspot", evidence2Source.sourceEvent());
        assertEquals(Sets.newHashSet("https://www.google.com/#q=FDA", "https://www.google.com/#q=NCCN"), evidence2Source.sourceUrls());
        assertEquals(ProtectEvidenceType.HOTSPOT_MUTATION, evidence2Source.evidenceType());
        assertNull(evidence2Source.rangeRank());
        assertTrue(evidence2Source.evidenceUrls().isEmpty());

        // evidences 3
        ProtectEvidence evidence3 = findByTreatmentAndEvent(evidences, "Dabrafenib + Trametinib", "p.Val600Glu");
        assertEquals(2, evidence3.sources().size());

        ProtectSource evidence3Source1 = findBySource(evidence3.sources(), Knowledgebase.VICC_CGI);
        assertEquals("hotspot", evidence3Source1.sourceEvent());
        assertTrue(evidence3Source1.sourceUrls().isEmpty());
        assertEquals(ProtectEvidenceType.HOTSPOT_MUTATION, evidence3Source1.evidenceType());
        assertNull(evidence3Source1.rangeRank());
        assertEquals(Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/25399551", "http://www.ncbi.nlm.nih.gov/pubmed/27283860"),
                evidence3Source1.evidenceUrls());

        ProtectSource evidence3Source2 = findBySource(evidence3.sources(), Knowledgebase.VICC_CIVIC);
        assertEquals("hotspot", evidence3Source2.sourceEvent());
        assertEquals(Sets.newHashSet("https://www.google.com/#q=FDA"), evidence3Source2.sourceUrls());
        assertEquals(ProtectEvidenceType.HOTSPOT_MUTATION, evidence3Source2.evidenceType());
        assertNull(evidence3Source2.rangeRank());
        assertTrue(evidence3Source2.evidenceUrls().isEmpty());

        // evidences 4
        ProtectEvidence evidence4 = findByTreatmentAndEvent(evidences, "Vemurafenib", "p.Val600Glu");
        assertEquals(1, evidence4.sources().size());

        ProtectSource evidence4Source = findBySource(evidence4.sources(), Knowledgebase.VICC_CGI);
        assertEquals("hotspot", evidence4Source.sourceEvent());
        assertEquals(Sets.newHashSet("https://www.google.com/#q=FDA"), evidence4Source.sourceUrls());
        assertEquals(ProtectEvidenceType.HOTSPOT_MUTATION, evidence4Source.evidenceType());
        assertEquals("600", String.valueOf(evidence4Source.rangeRank()));
        assertTrue(evidence4Source.evidenceUrls().isEmpty());

        // evidences 5
        ProtectEvidence evidence5 = findByTreatmentAndEvent(evidences, "Vemurafenib", "some mutation");
        assertEquals(1, evidence5.sources().size());

        ProtectSource evidence5Source1 = findBySource(evidence5.sources(), Knowledgebase.VICC_CGI);
        assertEquals("hotspot", evidence5Source1.sourceEvent());
        assertEquals(Sets.newHashSet("https://www.google.com/#q=FDA"), evidence5Source1.sourceUrls());
        assertEquals(ProtectEvidenceType.SIGNATURE, evidence5Source1.evidenceType());
        assertNull(evidence5Source1.rangeRank());
        assertTrue(evidence5Source1.evidenceUrls().isEmpty());
    }

    @Test
    public void canConvertSourcesBackAndForth() {
        Set<ProtectSource> sources = Sets.newHashSet();

        sources.add(ImmutableProtectSource.builder()
                .name(Knowledgebase.VICC_CGI)
                .sourceEvent("event 1")
                .sourceUrls(Sets.newHashSet("url1", "url2", "url3"))
                .evidenceType(ProtectEvidenceType.ANY_MUTATION)
                .rangeRank(1)
                .evidenceUrls(Sets.newHashSet("url4", "url5", "url6"))
                .build());

        sources.add(ImmutableProtectSource.builder()
                .name(Knowledgebase.VICC_CIVIC)
                .sourceEvent("event 2")
                .evidenceType(ProtectEvidenceType.HOTSPOT_MUTATION)
                .build());

        sources.add(ImmutableProtectSource.builder()
                .name(Knowledgebase.VICC_JAX)
                .sourceEvent("event 3")
                .evidenceType(ProtectEvidenceType.HOTSPOT_MUTATION)
                .evidenceUrls(Sets.newHashSet("url1"))
                .build());

        assertEquals(sources, ProtectEvidenceFile.stringToSources(ProtectEvidenceFile.sourcesToString(sources)));
    }

    @NotNull
    private static ProtectEvidence findByTreatmentAndEvent(@NotNull Iterable<ProtectEvidence> evidences, @NotNull String treatment,
            @NotNull String event) {
        for (ProtectEvidence evidence : evidences) {
            if (evidence.treatment().equals(treatment) && evidence.event().equals(event)) {
                return evidence;
            }
        }

        throw new IllegalStateException("Could not find evidence with treatment: " + treatment + " and event " + event);
    }

    @NotNull
    private static ProtectSource findBySource(@NotNull Iterable<ProtectSource> sources, @NotNull Knowledgebase sourceToFind) {
        for (ProtectSource source : sources) {
            if (source.name() == sourceToFind) {
                return source;
            }
        }

        throw new IllegalStateException("Could not find source: " + sourceToFind);
    }
}