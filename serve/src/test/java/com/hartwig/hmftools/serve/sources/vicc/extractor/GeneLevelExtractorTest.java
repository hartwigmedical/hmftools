package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.ONCO;
import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.TSG;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.actionability.gene.GeneLevelEvent;
import com.hartwig.hmftools.serve.gene.GeneLevelAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.check.GeneChecker;
import com.hartwig.hmftools.serve.sources.vicc.check.GeneCheckerTestFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GeneLevelExtractorTest {

    private static final GeneChecker HG19_GENE_CHECKER = GeneCheckerTestFactory.buildForHG19();

    @Test
    public void canExtractGeneLevelEventONCO() {
        GeneLevelExtractor geneLevelExtractor = createWithDriverGenes(createDriverGenes("STK11", "MET", "KIT"));
        GeneLevelAnnotation geneLevelEvent = geneLevelExtractor.extract("KIT", EventType.GENE_LEVEL, "KIT  positive");

        assertNotNull(geneLevelEvent);
        assertEquals("KIT", geneLevelEvent.gene());
        assertEquals(GeneLevelEvent.ACTIVATION, geneLevelEvent.event());
    }

    @Test
    public void canExtractGeneLevelEventTSG() {
        GeneLevelExtractor geneLevelExtractor = createWithDriverGenes(createDriverGenes("STK11", "MET", "KIT"));
        GeneLevelAnnotation geneLevelEvent = geneLevelExtractor.extract("TP53", EventType.GENE_LEVEL, "TP53  negative");

        assertNotNull(geneLevelEvent);
        assertEquals("TP53", geneLevelEvent.gene());
        assertEquals(GeneLevelEvent.INACTIVATION, geneLevelEvent.event());
    }

    @Test
    public void canExtractGeneLevelEventGeneral() {
        GeneLevelExtractor geneLevelExtractor = createWithDriverGenes(createDriverGenes("STK11", "MET", "KIT"));
        GeneLevelAnnotation geneLevelEvent = geneLevelExtractor.extract("STK11", EventType.GENE_LEVEL, "Truncating Mutations");

        assertNotNull(geneLevelEvent);
        assertEquals("STK11", geneLevelEvent.gene());
        assertEquals(GeneLevelEvent.INACTIVATION, geneLevelEvent.event());
    }

    @Test
    public void canExtractGeneLevelEventFusion() {
        GeneLevelExtractor geneLevelExtractor = createWithDriverGenes(createDriverGenes("STK11", "MET", "KIT"));
        GeneLevelAnnotation geneLevelEvent = geneLevelExtractor.extract("NTRK3", EventType.PROMISCUOUS_FUSION, "NTRK3 fusion");

        assertNotNull(geneLevelEvent);
        assertEquals("NTRK3", geneLevelEvent.gene());
        assertEquals(GeneLevelEvent.FUSION, geneLevelEvent.event());
    }

    @Test
    public void canExtractGeneLevelEvent() {
        List<DriverGene> driverGenes = createDriverGenes("NOTCH1", "NF2", "KRAS");

        assertEquals(GeneLevelEvent.ACTIVATION, GeneLevelExtractor.extractGeneLevelEvent("KRAS", driverGenes, "KRAS oncogenic mutation"));
        assertEquals(GeneLevelEvent.INACTIVATION, GeneLevelExtractor.extractGeneLevelEvent("NOTCH1", driverGenes, "LOSS-OF-FUNCTION"));
        assertEquals(GeneLevelEvent.ACTIVATION, GeneLevelExtractor.extractGeneLevelEvent("NF2", driverGenes, "MUTATION"));
        assertEquals(GeneLevelEvent.INACTIVATION, GeneLevelExtractor.extractGeneLevelEvent("NOTCH1", driverGenes, "MUTATION"));
        assertEquals(GeneLevelEvent.INACTIVATION, GeneLevelExtractor.extractGeneLevelEvent("NOTCH1", driverGenes, "NOTCH1 "));
        assertEquals(GeneLevelEvent.ANY_MUTATION, GeneLevelExtractor.extractGeneLevelEvent("BRCA1", driverGenes, "BRCA1"));

        assertNull(GeneLevelExtractor.extractGeneLevelEvent("KRAS", driverGenes, "not a gene level event"));
    }

    @Test
    public void canDetermineGeneLevelFromDriverGenes() {
        List<DriverGene> driverGenes = createDriverGenes("STK11", "MET", "KIT");

        assertEquals(GeneLevelEvent.ACTIVATION, GeneLevelExtractor.determineGeneLevelEventFromDriverGenes("MET", driverGenes));
        assertEquals(GeneLevelEvent.INACTIVATION, GeneLevelExtractor.determineGeneLevelEventFromDriverGenes("STK11", driverGenes));
        assertEquals(GeneLevelEvent.ANY_MUTATION, GeneLevelExtractor.determineGeneLevelEventFromDriverGenes("MAP1K1", driverGenes));
    }

    @NotNull
    private static GeneLevelExtractor createWithDriverGenes(@NotNull List<DriverGene> driverGenes) {
        return new GeneLevelExtractor(HG19_GENE_CHECKER, HG19_GENE_CHECKER, driverGenes);
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