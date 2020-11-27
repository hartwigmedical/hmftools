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
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GeneLevelExtractorTest {


    @Test
    public void canExtractGeneLevelEventGene(){
        List<DriverGene> driverGenes = createDriverGenes("STK11", "MET", "KIT");

        Feature featureOnco = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.CIVIC,
                "ENST00000318493",
                "MET",
                "MUTATION",
                "description",
                "chromosome",
                "pos").features().get(0);
        assertEquals(GeneLevelEvent.ACTIVATION, GeneLevelExtractor.extractGeneLevelEventGene(featureOnco, driverGenes));

        Feature featureTsg = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.CGI, null,
                "STK11",
                "STK11  mut",
                "description",
                "chromosome",
                "pos").features().get(0);
        assertEquals(GeneLevelEvent.INACTIVATION, GeneLevelExtractor.extractGeneLevelEventGene(featureTsg, driverGenes));

        Feature featureUnknown= ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.CIVIC,
                null,
                "MAP1K1",
                "Truncating Mutations",
                "description",
                "chromosome",
                "pos").features().get(0);
        assertEquals(GeneLevelEvent.UNKNOWN, GeneLevelExtractor.extractGeneLevelEventGene(featureUnknown, driverGenes));


    }

    @Test
    public void canExtractGeneLevelEvent(){

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

    @Test
    public void canExtractGeneLevelEventONCO() {
        Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature = Maps.newHashMap();

        GeneLevelExtractor geneLevelExtractor = new GeneLevelExtractor(HmfGenePanelSupplier.allGenesMap37(), createDriverGenes());
        ViccEntry viccEntry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.ONCOKB,
                "any",
                "KIT",
                "KIT  positive",
                "description",
                "chromosome",
                "pos");

        assertEquals(GeneLevelExtractor.extractGeneLevelEvent(viccEntry.features().get(0), createDriverGenes()), GeneLevelEvent.ACTIVATION);

        geneLevelEventsPerFeature.put(viccEntry.features().get(0),
                ImmutableGeneLevelAnnotation.builder()
                        .gene("KIT")
                        .event(GeneLevelEvent.ACTIVATION)
                        .build());
        assertEquals(geneLevelEventsPerFeature, geneLevelExtractor.extractGeneLevelEvents(viccEntry));
    }

    @Test
    public void canExtractGeneLevelEventTSG() {
        Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature = Maps.newHashMap();

        GeneLevelExtractor geneLevelExtractor = new GeneLevelExtractor(HmfGenePanelSupplier.allGenesMap37(), createDriverGenes());
        ViccEntry viccEntry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.ONCOKB,
                "any",
                "TP53",
                "TP53  negative",
                "description",
                "chromosome",
                "pos");

        assertEquals(GeneLevelExtractor.extractGeneLevelEvent(viccEntry.features().get(0), createDriverGenes()), GeneLevelEvent.INACTIVATION);

        geneLevelEventsPerFeature.put(viccEntry.features().get(0),
                ImmutableGeneLevelAnnotation.builder()
                        .gene("TP53")
                        .event(GeneLevelEvent.INACTIVATION)
                        .build());
        assertEquals(geneLevelEventsPerFeature, geneLevelExtractor.extractGeneLevelEvents(viccEntry));
    }

    @Test
    public void canExtractGeneLevelEventGeneral() {
        Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature = Maps.newHashMap();

        GeneLevelExtractor geneLevelExtractor = new GeneLevelExtractor(HmfGenePanelSupplier.allGenesMap37(), createDriverGenes());
        ViccEntry viccEntry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.ONCOKB,
                "any",
                "SRSF2",
                "Truncating Mutations",
                "description",
                "chromosome",
                "pos");

        assertEquals(GeneLevelExtractor.extractGeneLevelEvent(viccEntry.features().get(0), createDriverGenes()), GeneLevelEvent.INACTIVATION);

        geneLevelEventsPerFeature.put(viccEntry.features().get(0),
                ImmutableGeneLevelAnnotation.builder()
                        .gene("SRSF2")
                        .event(GeneLevelEvent.INACTIVATION)
                        .build());
        assertEquals(geneLevelEventsPerFeature, geneLevelExtractor.extractGeneLevelEvents(viccEntry));
    }

    @Test
    public void canExtractGeneLevelEventFusion() {
        Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature = Maps.newHashMap();

        GeneLevelExtractor geneLevelExtractor = new GeneLevelExtractor(HmfGenePanelSupplier.allGenesMap37(), createDriverGenes());
        ViccEntry viccEntry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.ONCOKB,
                "any",
                "NTRK3",
                "NTRK3 fusion",
                "description",
                "chromosome",
                "pos");

        geneLevelEventsPerFeature.put(viccEntry.features().get(0),
                ImmutableGeneLevelAnnotation.builder()
                        .gene("NTRK3")
                        .event(GeneLevelEvent.FUSION)
                        .build());
        assertEquals(geneLevelEventsPerFeature, geneLevelExtractor.extractGeneLevelEvents(viccEntry));
    }

    @NotNull
    private static List<DriverGene> createDriverGenes() {
        //TODO: Determine real TSG and ONCO genes
        DriverGene driverGeneTsg = ImmutableDriverGene.builder()
                .gene("SRSF2")
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

        DriverGene driverGeneOnco = ImmutableDriverGene.builder()
                .gene("KIT")
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
        return Lists.newArrayList(driverGeneTsg, driverGeneOnco);
    }
}