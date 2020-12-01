package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.ONCO;
import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.TSG;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.serve.actionability.gene.GeneLevelEvent;
import com.hartwig.hmftools.serve.sources.vicc.ViccTestFactory;
import com.hartwig.hmftools.serve.sources.vicc.annotation.GeneLevelAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.check.GeneChecker;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GeneLevelExtractorTest {

    @Test
    public void canExtractGeneLevelEventONCO() {
        GeneLevelExtractor geneLevelExtractor =
                new GeneLevelExtractor(HmfGenePanelSupplier.allGenesMap37(), createDriverGenes("STK11", "MET", "KIT"), new GeneChecker());
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("KIT", "KIT  positive");

        Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature = geneLevelExtractor.extractGeneLevelEvents(viccEntry);

        assertEquals(1, geneLevelEventsPerFeature.size());
        assertEquals("KIT", geneLevelEventsPerFeature.get(viccEntry.features().get(0)).gene());
        assertEquals(GeneLevelEvent.ACTIVATION, geneLevelEventsPerFeature.get(viccEntry.features().get(0)).event());
    }

    @Test
    public void canExtractGeneLevelEventTSG() {
        GeneLevelExtractor geneLevelExtractor =
                new GeneLevelExtractor(HmfGenePanelSupplier.allGenesMap37(), createDriverGenes("STK11", "MET", "KIT"), new GeneChecker());
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("TP53", "TP53  negative");

        Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature = geneLevelExtractor.extractGeneLevelEvents(viccEntry);

        assertEquals(1, geneLevelEventsPerFeature.size());
        assertEquals("TP53", geneLevelEventsPerFeature.get(viccEntry.features().get(0)).gene());
        assertEquals(GeneLevelEvent.INACTIVATION, geneLevelEventsPerFeature.get(viccEntry.features().get(0)).event());
    }

    @Test
    public void canExtractGeneLevelEventGeneral() {
        GeneLevelExtractor geneLevelExtractor =
                new GeneLevelExtractor(HmfGenePanelSupplier.allGenesMap37(), createDriverGenes("STK11", "MET", "KIT"), new GeneChecker());
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("STK11", "Truncating Mutations");

        Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature = geneLevelExtractor.extractGeneLevelEvents(viccEntry);

        assertEquals(1, geneLevelEventsPerFeature.size());
        assertEquals("STK11", geneLevelEventsPerFeature.get(viccEntry.features().get(0)).gene());
        assertEquals(GeneLevelEvent.INACTIVATION, geneLevelEventsPerFeature.get(viccEntry.features().get(0)).event());
    }

    @Test
    public void canExtractGeneLevelEventFusion() {
        GeneLevelExtractor geneLevelExtractor =
                new GeneLevelExtractor(HmfGenePanelSupplier.allGenesMap37(), createDriverGenes("STK11", "MET", "KIT"), new GeneChecker());
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("NTRK3", "NTRK3 fusion");

        Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature = geneLevelExtractor.extractGeneLevelEvents(viccEntry);

        assertEquals(1, geneLevelEventsPerFeature.size());
        assertEquals("NTRK3", geneLevelEventsPerFeature.get(viccEntry.features().get(0)).gene());
        assertEquals(GeneLevelEvent.FUSION, geneLevelEventsPerFeature.get(viccEntry.features().get(0)).event());
    }

    @Test
    public void canExtractGeneLevelEvent() {
        List<DriverGene> driverGenes = createDriverGenes("NOTCH1", "NF2", "KRAS");

        Feature featureActivation = ViccTestFactory.testFeatureWithGeneAndName("KRAS", "KRAS oncogenic mutation");
        assertEquals(GeneLevelEvent.ACTIVATION, GeneLevelExtractor.extractGeneLevelEvent(featureActivation, driverGenes));

        Feature featureInactivation = ViccTestFactory.testFeatureWithGeneAndName("NOTCH1", "LOSS-OF-FUNCTION");
        assertEquals(GeneLevelEvent.INACTIVATION, GeneLevelExtractor.extractGeneLevelEvent(featureInactivation, driverGenes));

        Feature featureGeneralActivation = ViccTestFactory.testFeatureWithGeneAndName("NF2", "MUTATION");
        assertEquals(GeneLevelEvent.ACTIVATION, GeneLevelExtractor.extractGeneLevelEvent(featureGeneralActivation, driverGenes));

        Feature featureGeneralInactivation = ViccTestFactory.testFeatureWithGeneAndName("NOTCH1", "MUTATION");
        assertEquals(GeneLevelEvent.INACTIVATION, GeneLevelExtractor.extractGeneLevelEvent(featureGeneralInactivation, driverGenes));

        Feature featureGeneOnly = ViccTestFactory.testFeatureWithGeneAndName("NOTCH1", "NOTCH1 ");
        assertEquals(GeneLevelEvent.INACTIVATION, GeneLevelExtractor.extractGeneLevelEvent(featureGeneOnly, driverGenes));

        Feature featureNoDriverGene = ViccTestFactory.testFeatureWithGeneAndName("BRCA1", "BRCA1");
        assertEquals(GeneLevelEvent.ANY_MUTATION, GeneLevelExtractor.extractGeneLevelEvent(featureNoDriverGene, driverGenes));

        Feature featureUnknown = ViccTestFactory.testFeatureWithGeneAndName("KRAS", "not a gene level event");
        assertNull(GeneLevelExtractor.extractGeneLevelEvent(featureUnknown, driverGenes));
    }

    @Test
    public void canDetermineGeneLevelFromDriverGenes() {
        List<DriverGene> driverGenes = createDriverGenes("STK11", "MET", "KIT");

        assertEquals(GeneLevelEvent.ACTIVATION, GeneLevelExtractor.determineGeneLevelEventFromDriverGenes("MET", driverGenes));
        assertEquals(GeneLevelEvent.INACTIVATION, GeneLevelExtractor.determineGeneLevelEventFromDriverGenes("STK11", driverGenes));
        assertEquals(GeneLevelEvent.ANY_MUTATION, GeneLevelExtractor.determineGeneLevelEventFromDriverGenes("MAP1K1", driverGenes));
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