package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.ONCO;
import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.TSG;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.serve.actionability.range.MutationTypeFilter;
import com.hartwig.hmftools.serve.exon.ExonAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.ViccTestFactory;
import com.hartwig.hmftools.serve.sources.vicc.check.GeneChecker;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ExonExtractorTest {

    private static final Map<String, HmfTranscriptRegion> HG19_GENE_MAP = HmfGenePanelSupplier.allGenesMap37();
    private static final GeneChecker HG19_GENE_CHECKER = new GeneChecker(HG19_GENE_MAP.keySet());

    @Test
    public void canExtractRangesExonAndFusion() {
        ExonExtractor extractor = createWithDriverGenes(createDriverGenes("TP53", "EGFR", "KIT"));
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("KIT", "EXON 11 MUTATION");

        Map<Feature, List<ExonAnnotation>> exonsPerFeature = extractor.extract(viccEntry);

        assertEquals(1, exonsPerFeature.size());
        assertEquals(1, exonsPerFeature.values().size());
        ExonAnnotation annotation = exonsPerFeature.get(viccEntry.features().get(0)).get(0);

        assertEquals("4", annotation.chromosome());
        assertEquals(55593577, annotation.start());
        assertEquals(55593713, annotation.end());
        assertEquals("KIT", annotation.gene());
        assertEquals(MutationTypeFilter.MISSENSE_ANY, annotation.mutationType());
    }
    
    @Test
    public void canExtractRangesExonForward() {
        ExonExtractor extractor = createWithDriverGenes(createDriverGenes("TP53", "EGFR", "KIT"));
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("EGFR", "EXON 19 DELETION");

        Map<Feature, List<ExonAnnotation>> exonsPerFeature = extractor.extract(viccEntry);

        assertEquals(1, exonsPerFeature.size());
        assertEquals(1, exonsPerFeature.values().size());
        ExonAnnotation annotation = exonsPerFeature.get(viccEntry.features().get(0)).get(0);

        assertEquals("7", annotation.chromosome());
        assertEquals(55242410, annotation.start());
        assertEquals(55242518, annotation.end());
        assertEquals("EGFR", annotation.gene());
        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_DELETION, annotation.mutationType());
    }

    @Test
    public void canExtractRangesExonReverse() {
        ExonExtractor extractor = createWithDriverGenes(createDriverGenes("TP53", "EGFR", "KIT"));
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("KRAS", "EXON 2 DELETION");

        Map<Feature, List<ExonAnnotation>> exonsPerFeature = extractor.extract(viccEntry);

        assertEquals(1, exonsPerFeature.size());
        assertEquals(1, exonsPerFeature.values().size());
        ExonAnnotation annotation = exonsPerFeature.get(viccEntry.features().get(0)).get(0);

        assertEquals("12", annotation.chromosome());
        assertEquals(25398203, annotation.start());
        assertEquals(25398334, annotation.end());
        assertEquals("KRAS", annotation.gene());
        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_DELETION, annotation.mutationType());
    }

    @Test
    public void canExtractAnnotationForEntryWithTwoFeaturesExon() {
        ExonExtractor extractor = createWithDriverGenes(Lists.newArrayList());

        Feature feature1 = ViccTestFactory.testFeatureWithGeneAndName("ERBB2", "Exon 20 insertions/deletions");
        Feature feature2 = ViccTestFactory.testFeatureWithGeneAndName("ERBB2", "Exon 20 insertions");
        ViccEntry entry = ViccTestFactory.testEntryWithFeatures(Lists.newArrayList(feature1, feature2));

        assertEquals(2, extractor.extract(entry).size());
        assertEquals(1, extractor.extract(entry).get(feature1).size());
        assertEquals(1, extractor.extract(entry).get(feature2).size());
    }

    @Test
    public void canExtractExons() {
        assertEquals(Lists.newArrayList(19), ExonExtractor.extractExonNumbers("EGFR exon 19 insertions"));
        assertEquals(Lists.newArrayList(20), ExonExtractor.extractExonNumbers("ERBB2 proximal exon 20"));
        assertEquals(Lists.newArrayList(9, 11, 13, 14, 17), ExonExtractor.extractExonNumbers("KIT mutation in exon 9,11,13,14 or 17"));
        assertEquals(Lists.newArrayList(16, 17, 18, 19), ExonExtractor.extractExonNumbers("MET mutation in exon 16-19"));
        assertEquals(Lists.newArrayList(2, 3), ExonExtractor.extractExonNumbers("Null (Partial deletion of Exons 2 & 3)"));
        assertEquals(Lists.newArrayList(12), ExonExtractor.extractExonNumbers("Exon 12 splice site insertion"));

        assertTrue(ExonExtractor.extractExonNumbers("Not an exon number").isEmpty());
    }

    @NotNull
    private static ExonExtractor createWithDriverGenes(@NotNull List<DriverGene> driverGenes) {
        return new ExonExtractor(HG19_GENE_CHECKER, driverGenes, HG19_GENE_MAP);
    }

    @NotNull
    private static List<DriverGene> createDriverGenes(@NotNull String geneTsg, @NotNull String geneOnco1, @NotNull String geneOnco2) {
        ImmutableDriverGene.Builder driverGeneBuilder = ImmutableDriverGene.builder()
                .reportMissenseAndInframe(false)
                .reportNonsenseAndFrameshift(false)
                .reportSplice(false)
                .reportDeletion(false)
                .reportDisruption(false)
                .reportAmplification(false)
                .reportSomaticHotspot(false)
                .reportGermlineVariant(false)
                .reportGermlineHotspot(false);

        DriverGene driverGeneTsg = driverGeneBuilder.gene(geneTsg).likelihoodType(TSG).build();
        DriverGene driverGeneOnco1 = driverGeneBuilder.gene(geneOnco1).likelihoodType(ONCO).build();
        DriverGene driverGeneOnco2 = driverGeneBuilder.gene(geneOnco2).likelihoodType(ONCO).build();

        return Lists.newArrayList(driverGeneTsg, driverGeneOnco1, driverGeneOnco2);
    }
}