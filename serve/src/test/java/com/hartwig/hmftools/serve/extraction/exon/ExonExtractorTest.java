package com.hartwig.hmftools.serve.extraction.exon;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.DriverGeneTestFactory;
import com.hartwig.hmftools.serve.EnsemblDataCacheTestFactory;
import com.hartwig.hmftools.serve.extraction.util.GeneChecker;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilter;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilterAlgo;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ExonExtractorTest {

    @Test
    public void canExtractExonForExonAndFusion() {
        ExonExtractor extractor = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("TP53", "KIT"));
        List<ExonAnnotation> exons = extractor.extract("KIT", null, EventType.FUSION_PAIR_AND_EXON, "EXON 11 MUTATION");

        assertEquals(1, exons.size());

        assertEquals("4", exons.get(0).chromosome());
        assertEquals(55593572, exons.get(0).start());
        assertEquals(55593718, exons.get(0).end());
        assertEquals("KIT", exons.get(0).gene());
        assertEquals(MutationTypeFilter.MISSENSE, exons.get(0).mutationType());
    }

    @Test
    public void canExtractExonForwardStrand() {
        ExonExtractor extractor = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("TP53", "EGFR"));
        List<ExonAnnotation> exons = extractor.extract("EGFR", null, EventType.EXON, "EXON 19 DELETION");

        assertEquals(1, exons.size());

        assertEquals("7", exons.get(0).chromosome());
        assertEquals(55242405, exons.get(0).start());
        assertEquals(55242523, exons.get(0).end());
        assertEquals("EGFR", exons.get(0).gene());
        assertEquals(MutationTypeFilter.INFRAME_DELETION, exons.get(0).mutationType());
    }

    @Test
    public void canExtractExonReverseStrand() {
        ExonExtractor extractor = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("TP53", "EGFR"));
        List<ExonAnnotation> exons = extractor.extract("KRAS", null, EventType.EXON, "EXON 2 DELETION");

        assertEquals(1, exons.size());

        assertEquals("12", exons.get(0).chromosome());
        assertEquals(25398198, exons.get(0).start());
        assertEquals(25398339, exons.get(0).end());
        assertEquals("KRAS", exons.get(0).gene());
        assertEquals(MutationTypeFilter.INFRAME_DELETION, exons.get(0).mutationType());
    }

    @Test
    public void canFilterOnNonCanonicalTranscript() {
        ExonExtractor extractor = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("TP53", "EGFR"));
        assertNull(extractor.extract("KRAS", "not the canonical transcript", EventType.EXON, "EXON 2 DELETION"));
    }

    @Test
    public void canFilterWhenExonIndicesDoNotExist() {
        ExonExtractor extractor = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("TP53", "EGFR"));
        assertNull(extractor.extract("KRAS", "ENST00000256078", EventType.EXON, "not a correct event"));
    }

    @Test
    public void canFilterWhenExonIndexNotOnTranscript() {
        ExonExtractor extractor = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("TP53", "EGFR"));
        assertNull(extractor.extract("KRAS", "ENST00000256078", EventType.EXON, "Exon 2000 deletion"));
    }

    @Test
    public void canExtractExonIndices() {
        assertEquals(Lists.newArrayList(19), ExonExtractor.extractExonIndices("EGFR exon 19 insertions"));
        assertEquals(Lists.newArrayList(20), ExonExtractor.extractExonIndices("ERBB2 proximal exon 20"));
        assertEquals(Lists.newArrayList(9, 11, 13, 14, 17), ExonExtractor.extractExonIndices("KIT mutation in exon 9,11,13,14 or 17"));
        assertEquals(Lists.newArrayList(16, 17, 18, 19), ExonExtractor.extractExonIndices("MET mutation in exon 16-19"));
        assertEquals(Lists.newArrayList(2, 3), ExonExtractor.extractExonIndices("Null (Partial deletion of Exons 2 & 3)"));
        assertEquals(Lists.newArrayList(12), ExonExtractor.extractExonIndices("Exon 12 splice site insertion"));

        assertNull(ExonExtractor.extractExonIndices("Not an exon number"));
    }

    @NotNull
    private static ExonExtractor createWithDriverGenes(@NotNull List<DriverGene> driverGenes) {
        return new ExonExtractor(new GeneChecker(Sets.newHashSet("TP53", "KIT", "EGFR", "KRAS")),
                new MutationTypeFilterAlgo(driverGenes),
                EnsemblDataCacheTestFactory.create37());
    }
}