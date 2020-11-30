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
import com.hartwig.hmftools.serve.actionability.gene.GeneLevelEvent;
import com.hartwig.hmftools.serve.sources.vicc.ViccTestFactory;
import com.hartwig.hmftools.serve.sources.vicc.annotation.GeneLevelAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.annotation.ImmutableGeneLevelAnnotation;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GeneLevelExtractorTest {

    @Test
    public void canExtractGeneLevelEventGene() {
        List<DriverGene> driverGenes = createDriverGenes("STK11", "MET", "KIT");

        Feature featureOnco = ViccTestFactory.testEntryWithGeneAndEvent("MET", "MUTATION").features().get(0);
        assertEquals(GeneLevelEvent.ACTIVATION, GeneLevelExtractor.extractGeneLevelEventGene(featureOnco, driverGenes));

        Feature featureTsg = ViccTestFactory.testEntryWithGeneAndEvent("STK11", "STK11  mut").features().get(0);
        assertEquals(GeneLevelEvent.INACTIVATION, GeneLevelExtractor.extractGeneLevelEventGene(featureTsg, driverGenes));

        Feature featureUnknown = ViccTestFactory.testEntryWithGeneAndEvent("MAP1K1", "Truncating Mutations").features().get(0);
        assertEquals(GeneLevelEvent.UNKNOWN, GeneLevelExtractor.extractGeneLevelEventGene(featureUnknown, driverGenes));
    }

    @Test
    public void canExtractGeneLevelEvent() {
        List<DriverGene> driverGenes = createDriverGenes("NOTCH1", "NF2", "KRAS");

        Feature featureActivation = ViccTestFactory.testEntryWithGeneAndEvent("KRAS", "KRAS oncogenic mutation").features().get(0);
        assertEquals(GeneLevelEvent.ACTIVATION, GeneLevelExtractor.extractGeneLevelEvent(featureActivation, driverGenes));

        Feature featureInactivation = ViccTestFactory.testEntryWithGeneAndEvent("NOTCH1", "LOSS-OF-FUNCTION").features().get(0);
        assertEquals(GeneLevelEvent.INACTIVATION, GeneLevelExtractor.extractGeneLevelEvent(featureInactivation, driverGenes));

        Feature featureGeneralActivation = ViccTestFactory.testEntryWithGeneAndEvent("NF2", "MUTATION").features().get(0);
        assertEquals(GeneLevelEvent.ACTIVATION, GeneLevelExtractor.extractGeneLevelEvent(featureGeneralActivation, driverGenes));

        Feature featureGeneralInactivation = ViccTestFactory.testEntryWithGeneAndEvent("NOTCH1", "MUTATION").features().get(0);
        assertEquals(GeneLevelEvent.INACTIVATION, GeneLevelExtractor.extractGeneLevelEvent(featureGeneralInactivation, driverGenes));

        Feature featureGeneOnly = ViccTestFactory.testEntryWithGeneAndEvent("NOTCH1", "NOTCH1 ").features().get(0);
        assertEquals(GeneLevelEvent.INACTIVATION, GeneLevelExtractor.extractGeneLevelEvent(featureGeneOnly, driverGenes));

        Feature featureUnkown = ViccTestFactory.testEntryWithGeneAndEvent("KRAS", "abcd").features().get(0);
        assertEquals(GeneLevelEvent.UNKNOWN, GeneLevelExtractor.extractGeneLevelEvent(featureUnkown, driverGenes));

    }

    @Test
    public void canExtractGeneLevelEventONCO() {
        Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature = Maps.newHashMap();

        GeneLevelExtractor geneLevelExtractor =
                new GeneLevelExtractor(HmfGenePanelSupplier.allGenesMap37(), createDriverGenes("STK11", "MET", "KIT"));
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("KIT", "KIT  positive");

        assertEquals(GeneLevelExtractor.extractGeneLevelEvent(viccEntry.features().get(0), createDriverGenes("STK11", "MET", "KIT")),
                GeneLevelEvent.ACTIVATION);

        geneLevelEventsPerFeature.put(viccEntry.features().get(0),
                ImmutableGeneLevelAnnotation.builder().gene("KIT").event(GeneLevelEvent.ACTIVATION).build());
        assertEquals(geneLevelEventsPerFeature, geneLevelExtractor.extractGeneLevelEvents(viccEntry));
    }

    @Test
    public void canExtractGeneLevelEventTSG() {
        Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature = Maps.newHashMap();

        GeneLevelExtractor geneLevelExtractor =
                new GeneLevelExtractor(HmfGenePanelSupplier.allGenesMap37(), createDriverGenes("STK11", "MET", "KIT"));
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("TP53", "TP53  negative");

        assertEquals(GeneLevelExtractor.extractGeneLevelEvent(viccEntry.features().get(0), createDriverGenes("STK11", "MET", "KIT")),
                GeneLevelEvent.INACTIVATION);

        geneLevelEventsPerFeature.put(viccEntry.features().get(0),
                ImmutableGeneLevelAnnotation.builder().gene("TP53").event(GeneLevelEvent.INACTIVATION).build());
        assertEquals(geneLevelEventsPerFeature, geneLevelExtractor.extractGeneLevelEvents(viccEntry));
    }

    @Test
    public void canExtractGeneLevelEventGeneral() {
        Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature = Maps.newHashMap();

        GeneLevelExtractor geneLevelExtractor =
                new GeneLevelExtractor(HmfGenePanelSupplier.allGenesMap37(), createDriverGenes("STK11", "MET", "KIT"));
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("STK11", "Truncating Mutations");

        assertEquals(GeneLevelExtractor.extractGeneLevelEvent(viccEntry.features().get(0), createDriverGenes("STK11", "MET", "KIT")),
                GeneLevelEvent.INACTIVATION);

        geneLevelEventsPerFeature.put(viccEntry.features().get(0),
                ImmutableGeneLevelAnnotation.builder().gene("STK11").event(GeneLevelEvent.INACTIVATION).build());
        assertEquals(geneLevelEventsPerFeature, geneLevelExtractor.extractGeneLevelEvents(viccEntry));
    }

    @Test
    public void canExtractGeneLevelEventFusion() {
        Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature = Maps.newHashMap();

        GeneLevelExtractor geneLevelExtractor =
                new GeneLevelExtractor(HmfGenePanelSupplier.allGenesMap37(), createDriverGenes("STK11", "MET", "KIT"));
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("NTRK3", "NTRK3 fusion");

        geneLevelEventsPerFeature.put(viccEntry.features().get(0),
                ImmutableGeneLevelAnnotation.builder().gene("NTRK3").event(GeneLevelEvent.FUSION).build());
        assertEquals(geneLevelEventsPerFeature, geneLevelExtractor.extractGeneLevelEvents(viccEntry));
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