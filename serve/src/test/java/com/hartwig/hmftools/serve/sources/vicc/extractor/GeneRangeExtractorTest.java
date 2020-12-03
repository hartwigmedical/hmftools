package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.ONCO;
import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.TSG;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.serve.actionability.range.MutationTypeFilter;
import com.hartwig.hmftools.serve.sources.vicc.ViccTestFactory;
import com.hartwig.hmftools.serve.sources.vicc.annotation.GeneRangeAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.annotation.GeneRangeType;
import com.hartwig.hmftools.serve.sources.vicc.annotation.ImmutableGeneRangeAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.check.GeneChecker;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GeneRangeExtractorTest {

    private static final Map<String, HmfTranscriptRegion> HG19_GENE_MAP = HmfGenePanelSupplier.allGenesMap37();
    private static final GeneChecker HG19_GENE_CHECKER = new GeneChecker(HG19_GENE_MAP.keySet());

    @Test
    public void canExtractRangesExonAndFusion() {
        GeneRangeExtractor geneRangeExtractor = createWithDriverGenes(createDriverGenes("TP53", "EGFR", "KIT"));
        Map<Feature, List<GeneRangeAnnotation>> geneRangesPerFeature = Maps.newHashMap();
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("KIT", "EXON 11 MUTATION");

        geneRangesPerFeature.put(viccEntry.features().get(0),
                Lists.newArrayList(ImmutableGeneRangeAnnotation.builder()
                        .gene("KIT")
                        .chromosome("4")
                        .start(55593577)
                        .end(55593713)
                        .mutationType(MutationTypeFilter.MISSENSE_ANY)
                        .rangeType(GeneRangeType.EXON)
                        .rangeNumber(11)
                        .exonId("ENSE00001074417")
                        .geneOrientation(Strand.FORWARD)
                        .build()));

        assertEquals(geneRangesPerFeature, geneRangeExtractor.extractGeneRanges(viccEntry));
    }

    @Test
    public void canExtractRangesCodonOnMultipleExons() {
        GeneRangeExtractor geneRangeExtractor = createWithDriverGenes(createDriverGenes("TP53", "KRAS", "KIT"));
        Map<Feature, List<GeneRangeAnnotation>> geneRangesPerFeature = Maps.newHashMap();
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("KRAS", "R97");

        geneRangesPerFeature.put(viccEntry.features().get(0),
                Lists.newArrayList(ImmutableGeneRangeAnnotation.builder()
                        .gene("KRAS")
                        .chromosome("12")
                        .start(25380168)
                        .end(25378707)
                        .mutationType(MutationTypeFilter.MISSENSE_ANY)
                        .rangeType(GeneRangeType.CODON)
                        .rangeNumber(97)
                        .geneOrientation(Strand.REVERSE)
                        .build()));

        assertEquals(geneRangesPerFeature, geneRangeExtractor.extractGeneRanges(viccEntry));
    }

    @Test
    public void canExtractRangesCodon() {
        GeneRangeExtractor geneRangeExtractor = createWithDriverGenes(createDriverGenes("TP53", "EGFR", "KIT"));
        Map<Feature, List<GeneRangeAnnotation>> geneRangesPerFeature = Maps.newHashMap();
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("TP53", "R249");

        geneRangesPerFeature.put(viccEntry.features().get(0),
                Lists.newArrayList(ImmutableGeneRangeAnnotation.builder()
                        .gene("TP53")
                        .chromosome("17")
                        .start(7577536)
                        .end(7577534)
                        .mutationType(MutationTypeFilter.ANY)
                        .rangeType(GeneRangeType.CODON)
                        .rangeNumber(249)
                        .geneOrientation(Strand.REVERSE)
                        .build()));

        assertEquals(geneRangesPerFeature, geneRangeExtractor.extractGeneRanges(viccEntry));
    }

    @Test
    public void canExtractRangesExonForward() {
        GeneRangeExtractor geneRangeExtractor = createWithDriverGenes(createDriverGenes("TP53", "EGFR", "KIT"));
        Map<Feature, List<GeneRangeAnnotation>> geneRangesPerFeature = Maps.newHashMap();
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("EGFR", "EXON 19 DELETION");

        geneRangesPerFeature.put(viccEntry.features().get(0),
                Lists.newArrayList(ImmutableGeneRangeAnnotation.builder()
                        .gene("EGFR")
                        .chromosome("7")
                        .start(55242410)
                        .end(55242518)
                        .mutationType(MutationTypeFilter.MISSENSE_INFRAME_DELETION)
                        .rangeType(GeneRangeType.EXON)
                        .rangeNumber(19)
                        .exonId("ENSE00001756460")
                        .geneOrientation(Strand.FORWARD)
                        .build()));

        assertEquals(geneRangesPerFeature, geneRangeExtractor.extractGeneRanges(viccEntry));
    }

    @Test
    public void canExtractRangesExonReverse() {
        GeneRangeExtractor geneRangeExtractor = createWithDriverGenes(createDriverGenes("TP53", "EGFR", "KIT"));
        Map<Feature, List<GeneRangeAnnotation>> geneRangesPerFeature = Maps.newHashMap();
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("KRAS", "EXON 2 DELETION");

        geneRangesPerFeature.put(viccEntry.features().get(0),
                Lists.newArrayList(ImmutableGeneRangeAnnotation.builder()
                        .gene("KRAS")
                        .chromosome("12")
                        .start(25398203)
                        .end(25398334)
                        .mutationType(MutationTypeFilter.MISSENSE_INFRAME_DELETION)
                        .rangeType(GeneRangeType.EXON)
                        .rangeNumber(2)
                        .exonId("ENSE00000936617")
                        .geneOrientation(Strand.REVERSE)
                        .build()));

        assertEquals(geneRangesPerFeature, geneRangeExtractor.extractGeneRanges(viccEntry));
    }

    @Test
    public void canExtractAnnotationForEntryWithTwoFeaturesExon() {
        GeneRangeExtractor extractor = createWithDriverGenes(Lists.newArrayList());

        Feature feature1 = ViccTestFactory.testFeatureWithGeneAndName("ERBB2", "Exon 20 insertions/deletions");
        Feature feature2 = ViccTestFactory.testFeatureWithGeneAndName("ERBB2", "Exon 20 insertions");
        ViccEntry entry = ViccTestFactory.testEntryWithFeatures(Lists.newArrayList(feature1, feature2));

        assertEquals(2, extractor.extractGeneRanges(entry).size());
        assertEquals(1, extractor.extractGeneRanges(entry).get(feature1).size());
        assertEquals(1, extractor.extractGeneRanges(entry).get(feature2).size());
    }

    @Test
    public void canExtractAnnotationForEntryWithTwoFeaturesCodon() {
        GeneRangeExtractor extractor = createWithDriverGenes(createDriverGenes("TP53", "PIK3CA", "KRAS"));

        Feature feature1 = ViccTestFactory.testFeatureWithGeneAndName("PIK3CA", "E545X");
        Feature feature2 = ViccTestFactory.testFeatureWithGeneAndName("KRAS", "G12X");
        ViccEntry entry = ViccTestFactory.testEntryWithFeatures(Lists.newArrayList(feature1, feature2));

        assertEquals(2, extractor.extractGeneRanges(entry).size());
        assertEquals(1, extractor.extractGeneRanges(entry).get(feature1).size());
        assertEquals(1, extractor.extractGeneRanges(entry).get(feature2).size());
    }

    @Test
    public void canExtractCodons() {
        assertEquals(Integer.valueOf("600"), GeneRangeExtractor.extractCodonNumber("BRAF (V600)"));
        assertEquals(Integer.valueOf("742"), GeneRangeExtractor.extractCodonNumber("W742"));
        assertEquals(Integer.valueOf("179"), GeneRangeExtractor.extractCodonNumber("Q179X"));
        assertEquals(Integer.valueOf("61"), GeneRangeExtractor.extractCodonNumber("KRAS Q61X"));

        assertNull(GeneRangeExtractor.extractCodonNumber("Not a codon number"));
    }

    @Test
    public void canExtractExons() {
        assertEquals(Lists.newArrayList(19), GeneRangeExtractor.extractExonNumbers("EGFR exon 19 insertions"));
        assertEquals(Lists.newArrayList(20), GeneRangeExtractor.extractExonNumbers("ERBB2 proximal exon 20"));
        assertEquals(Lists.newArrayList(9, 11, 13, 14, 17), GeneRangeExtractor.extractExonNumbers("KIT mutation in exon 9,11,13,14 or 17"));
        assertEquals(Lists.newArrayList(16, 17, 18, 19), GeneRangeExtractor.extractExonNumbers("MET mutation in exon 16-19"));
        assertEquals(Lists.newArrayList(2, 3), GeneRangeExtractor.extractExonNumbers("Null (Partial deletion of Exons 2 & 3)"));
        assertEquals(Lists.newArrayList(12), GeneRangeExtractor.extractExonNumbers("Exon 12 splice site insertion"));

        assertTrue(GeneRangeExtractor.extractExonNumbers("Not an exon number").isEmpty());
    }

    @Test
    public void canExtractMutationFilter() {
        List<DriverGene> driverGenes = createDriverGenes("TP53", "EGFR", "ERBB2");
        String gene = "ERBB2";

        assertEquals(MutationTypeFilter.NONSENSE_OR_FRAMESHIFT,
                GeneRangeExtractor.extractMutationTypeFilter("EXON 9 FRAMESHIFT", driverGenes, gene));

        assertEquals(MutationTypeFilter.SPLICE,
                GeneRangeExtractor.extractMutationTypeFilter("Exon 12 splice site insertion", driverGenes, gene));

        assertEquals(MutationTypeFilter.SPLICE,
                GeneRangeExtractor.extractMutationTypeFilter("EXON 14 SKIPPING MUTATION", driverGenes, gene));

        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_DELETION,
                GeneRangeExtractor.extractMutationTypeFilter("EGFR exon 19 deletions", driverGenes, gene));

        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_DELETION,
                GeneRangeExtractor.extractMutationTypeFilter("Null (Partial deletion of Exons 2 & 3)", driverGenes, gene));

        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_INSERTION,
                GeneRangeExtractor.extractMutationTypeFilter("Exon 20 insertions", driverGenes, gene));

        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_ANY,
                GeneRangeExtractor.extractMutationTypeFilter("Exon 20 insertions/deletions", driverGenes, gene));

        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_ANY,
                GeneRangeExtractor.extractMutationTypeFilter("Exon 19 deletion/insertion", driverGenes, gene));

        assertEquals(MutationTypeFilter.UNKNOWN, GeneRangeExtractor.extractMutationTypeFilter("abcd", driverGenes, "efgh"));

        assertEquals(MutationTypeFilter.MISSENSE_ANY, GeneRangeExtractor.extractMutationTypeFilter("mut", driverGenes, gene));

        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_INSERTION,
                GeneRangeExtractor.extractMutationTypeFilter("insertion", driverGenes, gene));

        assertEquals(MutationTypeFilter.ANY, GeneRangeExtractor.extractMutationTypeFilter("mut", driverGenes, "TP53"));

        assertEquals(MutationTypeFilter.NONSENSE_OR_FRAMESHIFT,
                GeneRangeExtractor.extractMutationTypeFilter("frameshift", driverGenes, "TP53"));
    }

    @NotNull
    private static GeneRangeExtractor createWithDriverGenes(@NotNull List<DriverGene> driverGenes) {
        return new GeneRangeExtractor(HG19_GENE_CHECKER, driverGenes, HG19_GENE_MAP);
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