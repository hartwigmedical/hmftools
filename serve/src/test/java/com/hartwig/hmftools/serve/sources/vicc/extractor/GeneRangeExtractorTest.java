package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.ONCO;
import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.TSG;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.serve.actionability.range.MutationTypeFilter;
import com.hartwig.hmftools.serve.sources.vicc.ViccTestFactory;
import com.hartwig.hmftools.serve.sources.vicc.annotation.GeneRangeAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.annotation.ImmutableGeneRangeAnnotation;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GeneRangeExtractorTest {

    @Test
    public void canExtractCodon() {
        assertEquals(Integer.valueOf("600"), GeneRangeExtractor.extractCodonNumber("BRAF (V600)"));
        assertEquals(Integer.valueOf("742"), GeneRangeExtractor.extractCodonNumber("W742"));
        assertEquals(Integer.valueOf("179"), GeneRangeExtractor.extractCodonNumber("Q179X"));
        assertEquals(Integer.valueOf("61"), GeneRangeExtractor.extractCodonNumber("KRAS Q61X"));
    }

    @Test
    public void canExtractExon() {
        assertEquals(Lists.newArrayList(19), GeneRangeExtractor.extractExonNumbers("EGFR exon 19 insertions"));
        assertEquals(Lists.newArrayList(20), GeneRangeExtractor.extractExonNumbers("ERBB2 proximal exon 20"));
        assertEquals(Lists.newArrayList(9, 11, 13, 14, 17), GeneRangeExtractor.extractExonNumbers("KIT mutation in exon 9,11,13,14 or 17"));
        assertEquals(Lists.newArrayList(16, 17, 18, 19), GeneRangeExtractor.extractExonNumbers("MET mutation in exon 16-19"));
        assertEquals(Lists.newArrayList(2, 3), GeneRangeExtractor.extractExonNumbers("Null (Partial deletion of Exons 2 & 3)"));
        assertEquals(Lists.newArrayList(12), GeneRangeExtractor.extractExonNumbers("Exon 12 splice site insertion"));
    }

    @Test
    public void canExtractAnnotationForEntryWithTwoFeatures() {
        GeneRangeExtractor extractor = new GeneRangeExtractor(HmfGenePanelSupplier.allGenesMap37(), Lists.newArrayList());

        Feature feature1 = ViccTestFactory.testFeatureWithGeneAndName("ERBB2", "Exon 20 insertions/deletions");
        Feature feature2 = ViccTestFactory.testFeatureWithGeneAndName("ERBB2", "Exon 20 insertions");
        ViccEntry entry = ViccTestFactory.testEntryWithFeatures(Lists.newArrayList(feature1, feature2));

        assertEquals(2, extractor.extractGeneRanges(entry).size());
        assertEquals(1, extractor.extractGeneRanges(entry).get(feature1).size());
        assertEquals(1, extractor.extractGeneRanges(entry).get(feature2).size());
    }

    @Test
    public void canExtractMutationFilter() {
        Feature featureNonsenseFrameshift = ViccTestFactory.testEntryWithGeneAndEvent("CALR", "EXON 9 FRAMESHIFT").features().get(0);
        assertEquals(MutationTypeFilter.NONSENSE_OR_FRAMESHIFT,
                GeneRangeExtractor.extractSpecificMutationTypeFilter(featureNonsenseFrameshift));

        Feature featureSplice1 = ViccTestFactory.testEntryWithGeneAndEvent("JAK2", "Exon 12 splice site insertion").features().get(0);
        assertEquals(MutationTypeFilter.SPLICE, GeneRangeExtractor.extractSpecificMutationTypeFilter(featureSplice1));

        Feature featureSplice2 = ViccTestFactory.testEntryWithGeneAndEvent("MET", "EXON 14 SKIPPING MUTATION").features().get(0);
        assertEquals(MutationTypeFilter.SPLICE, GeneRangeExtractor.extractSpecificMutationTypeFilter(featureSplice2));

        Feature featureMissenseInframeDeletion1 =
                ViccTestFactory.testEntryWithGeneAndEvent("EGFR", "EGFR exon 19 deletions").features().get(0);
        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_DELETION,
                GeneRangeExtractor.extractSpecificMutationTypeFilter(featureMissenseInframeDeletion1));

        Feature featureMissenseInframeDeletion2 =
                ViccTestFactory.testEntryWithGeneAndEvent("VHL", "Null (Partial deletion of Exons 2 & 3)").features().get(0);
        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_DELETION,
                GeneRangeExtractor.extractSpecificMutationTypeFilter(featureMissenseInframeDeletion2));

        Feature featureMissenseInframeInsertion =
                ViccTestFactory.testEntryWithGeneAndEvent("ERBB2", "Exon 20 insertions").features().get(0);
        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_INSERTION,
                GeneRangeExtractor.extractSpecificMutationTypeFilter(featureMissenseInframeInsertion));

        Feature featureMissenseInframeAny1 =
                ViccTestFactory.testEntryWithGeneAndEvent("ERBB2", "Exon 20 insertions/deletions").features().get(0);
        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_ANY,
                GeneRangeExtractor.extractSpecificMutationTypeFilter(featureMissenseInframeAny1));

        Feature featureMissenseInframeAny2 =
                ViccTestFactory.testEntryWithGeneAndEvent("EGFR", "Exon 19 deletion/insertion").features().get(0);
        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_ANY,
                GeneRangeExtractor.extractSpecificMutationTypeFilter(featureMissenseInframeAny2));

        Feature featureUnknown = ViccTestFactory.testEntryWithGeneAndEvent("GNAQ", "abcd").features().get(0);
        assertEquals(MutationTypeFilter.UNKNOWN, GeneRangeExtractor.extractSpecificMutationTypeFilter(featureUnknown));
    }

    @Test
    public void canExtractSpecificMutationTypeFilter() {
        Feature featureMutationUnknownOnco = ViccTestFactory.testEntryWithGeneAndEvent("ERBB2", "D835").features().get(0);
        List<DriverGene> driverGenes = createDriverGenes("TP53", "EGFR", "ERBB2");
        String gene = "ERBB2";
        MutationTypeFilter filterUnknownOnco = MutationTypeFilter.UNKNOWN;
        assertEquals(MutationTypeFilter.MISSENSE_ANY,
                GeneRangeExtractor.extractMutationFilter(driverGenes, gene, filterUnknownOnco, featureMutationUnknownOnco));

        Feature featureMutationKnownOnco = ViccTestFactory.testEntryWithGeneAndEvent("ERBB2", "D835").features().get(0);
        MutationTypeFilter filterKnownOnco = MutationTypeFilter.MISSENSE_INFRAME_INSERTION;
        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_INSERTION,
                GeneRangeExtractor.extractMutationFilter(driverGenes, gene, filterKnownOnco, featureMutationKnownOnco));

        Feature featureMutationUnknownTsg = ViccTestFactory.testEntryWithGeneAndEvent("TP53", "R249").features().get(0);
        MutationTypeFilter filterUnknownTsg = MutationTypeFilter.UNKNOWN;
        assertEquals(MutationTypeFilter.MISSENSE_ANY,
                GeneRangeExtractor.extractMutationFilter(driverGenes, gene, filterUnknownTsg, featureMutationUnknownTsg));

        Feature featureMutationKnownTsg = ViccTestFactory.testEntryWithGeneAndEvent("TP53", "R249").features().get(0);
        MutationTypeFilter filterKnownTsg = MutationTypeFilter.NONSENSE_OR_FRAMESHIFT;
        assertEquals(MutationTypeFilter.NONSENSE_OR_FRAMESHIFT,
                GeneRangeExtractor.extractMutationFilter(driverGenes, gene, filterKnownTsg, featureMutationKnownTsg));
    }

    @Test
    public void canExtractRangesCodon() {
        GeneRangeExtractor geneRangeExtractor =
                new GeneRangeExtractor(HmfGenePanelSupplier.allGenesMap37(), createDriverGenes("TP53", "EGFR", "KIT"));
        Map<Feature, List<GeneRangeAnnotation>> geneRangesPerFeature = Maps.newHashMap();
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("TP53", "R249");

        geneRangesPerFeature.put(viccEntry.features().get(0),
                Lists.newArrayList(ImmutableGeneRangeAnnotation.builder()
                        .gene("TP53")
                        .chromosome("17")
                        .start(7577534)
                        .end(7577536)
                        .rangeInfo(249)
                        .mutationType(MutationTypeFilter.ANY)
                        .build()));

        assertEquals(geneRangesPerFeature, geneRangeExtractor.extractGeneRanges(viccEntry));
    }

    @Test
    public void canExtractRangesExon() {
        GeneRangeExtractor geneRangeExtractor =
                new GeneRangeExtractor(HmfGenePanelSupplier.allGenesMap37(), createDriverGenes("TP53", "EGFR", "KIT"));
        Map<Feature, List<GeneRangeAnnotation>> geneRangesPerFeature = Maps.newHashMap();
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("EGFR", "EXON 19 DELETION");

        geneRangesPerFeature.put(viccEntry.features().get(0),
                Lists.newArrayList(ImmutableGeneRangeAnnotation.builder()
                        .gene("EGFR")
                        .chromosome("7")
                        .start(55242410)
                        .end(55242518)
                        .rangeInfo(19)
                        .exonId("ENSE00001756460")
                        .mutationType(MutationTypeFilter.MISSENSE_INFRAME_DELETION)
                        .build()));

        assertEquals(geneRangesPerFeature, geneRangeExtractor.extractGeneRanges(viccEntry));
    }

    @Test
    public void canExtractRangesExonAndFusion() {
        GeneRangeExtractor geneRangeExtractor =
                new GeneRangeExtractor(HmfGenePanelSupplier.allGenesMap37(), createDriverGenes("TP53", "EGFR", "KIT"));
        Map<Feature, List<GeneRangeAnnotation>> geneRangesPerFeature = Maps.newHashMap();
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("KIT", "EXON 11 MUTATION");

        geneRangesPerFeature.put(viccEntry.features().get(0),
                Lists.newArrayList(ImmutableGeneRangeAnnotation.builder()
                        .gene("KIT")
                        .chromosome("4")
                        .start(55593577)
                        .end(55593713)
                        .rangeInfo(11)
                        .exonId("ENSE00001074417")
                        .mutationType(MutationTypeFilter.MISSENSE_ANY)
                        .build()));

        assertEquals(geneRangesPerFeature, geneRangeExtractor.extractGeneRanges(viccEntry));
    }

    @NotNull
    private static List<DriverGene> createDriverGenes(@NotNull String geneTsg, @NotNull String geneOnco1, @NotNull String geneOnco2) {
        DriverGene driverGeneTsg = ImmutableDriverGene.builder()
                .gene(geneTsg)
                .reportMissenseAndInframe(true)
                .reportNonsenseAndFrameshift(true)
                .reportSplice(true)
                .reportDeletion(true)
                .reportDisruption(true)
                .reportAmplification(true)
                .reportSomaticHotspot(true)
                .reportGermlineVariant(false)
                .reportGermlineHotspot(false)
                .likelihoodType(TSG)
                .build();

        DriverGene driverGeneOnco1 = ImmutableDriverGene.builder()
                .gene(geneOnco1)
                .reportMissenseAndInframe(true)
                .reportNonsenseAndFrameshift(true)
                .reportSplice(true)
                .reportDeletion(true)
                .reportDisruption(true)
                .reportAmplification(true)
                .reportSomaticHotspot(true)
                .reportGermlineVariant(false)
                .reportGermlineHotspot(false)
                .likelihoodType(ONCO)
                .build();

        DriverGene driverGeneOnco2 = ImmutableDriverGene.builder()
                .gene(geneOnco2)
                .reportMissenseAndInframe(true)
                .reportNonsenseAndFrameshift(true)
                .reportSplice(true)
                .reportDeletion(true)
                .reportDisruption(true)
                .reportAmplification(true)
                .reportSomaticHotspot(true)
                .reportGermlineVariant(false)
                .reportGermlineHotspot(false)
                .likelihoodType(ONCO)
                .build();
        return Lists.newArrayList(driverGeneTsg, driverGeneOnco1, driverGeneOnco2);
    }
}