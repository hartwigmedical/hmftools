package com.hartwig.hmftools.common.protect;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
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

        //evidences1
        ProtectEvidence evidence1 = findByTreatment(evidences, "Cobimetinib + Vemurafenib", "p.Val600Glu");
        assertEquals(1, evidence1.protectSources().size());

        List<ProtectSource> evidence1_source1List = Lists.newArrayList(evidence1.protectSources());
        ProtectSource evidence1_source1 = findBySource(evidence1_source1List, Knowledgebase.VICC_CGI);

        assertEquals(Knowledgebase.VICC_CGI, evidence1_source1.sources());
        assertEquals("hotspot", evidence1_source1.sourceEvent());
        assertEquals(Sets.newHashSet("https://www.google.com/#q=FDA"),
                evidence1_source1.sourceUrls());
        assertEquals(ProtectEvidenceType.HOTSPOT_MUTATION, evidence1_source1.evidenceType());
        assertNull(evidence1_source1.rangeRank());

        //evidences2
        ProtectEvidence evidence2 = findByTreatment(evidences, "Dabrafenib", "p.Val600Glu");
        assertEquals(1, evidence2.protectSources().size());

        List<ProtectSource> evidence2_source1List = Lists.newArrayList(evidence2.protectSources());
        ProtectSource evidence2_source1 = findBySource(evidence2_source1List, Knowledgebase.VICC_CGI);

        assertEquals(Knowledgebase.VICC_CGI, evidence2_source1.sources());
        assertEquals("hotspot", evidence2_source1.sourceEvent());
        assertEquals(Sets.newHashSet("https://www.google.com/#q=FDA", "https://www.google.com/#q=NCCN"),
                evidence2_source1.sourceUrls());
        assertEquals(ProtectEvidenceType.HOTSPOT_MUTATION, evidence2_source1.evidenceType());
        assertNull(evidence2_source1.rangeRank());

        //evidences3
        ProtectEvidence evidence3 = findByTreatment(evidences, "Dabrafenib + Trametinib", "p.Val600Glu");
        assertEquals(2, evidence3.protectSources().size());

        List<ProtectSource> evidence3_source1List = Lists.newArrayList(evidence3.protectSources());
        ProtectSource source3_evidence1 = findBySource(evidence3_source1List, Knowledgebase.VICC_CGI);

        assertEquals(Knowledgebase.VICC_CGI, source3_evidence1.sources());
        assertEquals("hotspot", source3_evidence1.sourceEvent());
        assertEquals(Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/25399551", "http://www.ncbi.nlm.nih.gov/pubmed/27283860"),
                source3_evidence1.sourceUrls());
        assertEquals(ProtectEvidenceType.HOTSPOT_MUTATION, source3_evidence1.evidenceType());
        assertNull(source3_evidence1.rangeRank());

        ProtectSource source3_evidence2 = findBySource(evidence3_source1List, Knowledgebase.VICC_CIVIC);
        assertEquals(Knowledgebase.VICC_CIVIC, source3_evidence2.sources());
        assertEquals("hotspot", source3_evidence2.sourceEvent());
        assertEquals(Sets.newHashSet("https://www.google.com/#q=FDA"),
                source3_evidence2.sourceUrls());
        assertEquals(ProtectEvidenceType.HOTSPOT_MUTATION, source3_evidence2.evidenceType());
        assertNull(source3_evidence2.rangeRank());

        //evidences4
        ProtectEvidence evidence4 = findByTreatment(evidences, "Vemurafenib", "p.Val600Glu");
        assertEquals(1, evidence4.protectSources().size());

        List<ProtectSource> evidence4_source_1_list = Lists.newArrayList(evidence4.protectSources());
        ProtectSource evidence4_source_1 = findBySource(evidence4_source_1_list, Knowledgebase.VICC_CGI);

        assertEquals(Knowledgebase.VICC_CGI, evidence4_source_1.sources());
        assertEquals("hotspot", evidence4_source_1.sourceEvent());
        assertEquals(Sets.newHashSet("https://www.google.com/#q=FDA"),
                evidence4_source_1.sourceUrls());
        assertEquals(ProtectEvidenceType.HOTSPOT_MUTATION, evidence4_source_1.evidenceType());
        assertEquals("600", String.valueOf(evidence4_source_1.rangeRank()));

        //evidences5
        ProtectEvidence evidence5 = findByTreatment(evidences, "Vemurafenib", "some mutation");
        assertEquals(1, evidence5.protectSources().size());

        List<ProtectSource> evidence5_source_1_list = Lists.newArrayList(evidence5.protectSources());
        ProtectSource evidence5_source1 = findBySource(evidence5_source_1_list, Knowledgebase.VICC_CGI);

        assertEquals(Knowledgebase.VICC_CGI, evidence5_source1.sources());
        assertEquals("hotspot", evidence5_source1.sourceEvent());
        assertEquals(Sets.newHashSet("https://www.google.com/#q=FDA"),
                evidence5_source1.sourceUrls());
        assertEquals(ProtectEvidenceType.SIGNATURE, evidence5_source1.evidenceType());
        assertNull(evidence5_source1.rangeRank());
    }

    @NotNull
    private static ProtectEvidence findByTreatment(@NotNull List<ProtectEvidence> evidences, @NotNull String treatment, @NotNull String event) {
        for (ProtectEvidence evidence : evidences) {
            if (evidence.treatment().equals(treatment) && evidence.event().equals(event)) {
                return evidence;
            }
        }

        throw new IllegalStateException("Could not find evidence with treatment: " + treatment);
    }

    @NotNull
    private static ProtectSource findBySource(@NotNull List<ProtectSource> sources, @NotNull Knowledgebase evidenceSource) {
        for (ProtectSource source : sources) {
            if (source.sources() == evidenceSource) {
                return source;
            }
        }

        throw new IllegalStateException("Could not find evidence of source: " + evidenceSource);
    }
}