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
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.serve.actionability.range.MutationTypeFilter;
import com.hartwig.hmftools.serve.codon.CodonAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.ViccTestFactory;
import com.hartwig.hmftools.serve.sources.vicc.check.GeneChecker;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CodonExtractorTest {

    private static final Map<String, HmfTranscriptRegion> HG19_GENE_MAP = HmfGenePanelSupplier.allGenesMap37();
    private static final GeneChecker HG19_GENE_CHECKER = new GeneChecker(HG19_GENE_MAP.keySet());

    @Test
    public void canExtractSimpleCodon() {
        CodonExtractor extractor = createWithDriverGenes(createDriverGenes("TP53", "EGFR", "KIT"));
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("TP53", "R249");

        Map<Feature, List<CodonAnnotation>> codonsPerFeature = extractor.extract(viccEntry);

        assertEquals(1, codonsPerFeature.size());
        assertEquals(1, codonsPerFeature.values().size());
        CodonAnnotation annotation = codonsPerFeature.get(viccEntry.features().get(0)).get(0);

        assertEquals("17", annotation.chromosome());
        assertEquals(7577534, annotation.start());
        assertEquals(7577536, annotation.end());
        assertEquals("TP53", annotation.gene());
        assertEquals(MutationTypeFilter.ANY, annotation.mutationType());
    }

    @Test
    public void canExtractCodonOnMultipleExons() {
        CodonExtractor extractor = createWithDriverGenes(createDriverGenes("TP53", "KRAS", "KIT"));
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("KRAS", "R97");

        Map<Feature, List<CodonAnnotation>> codonsPerFeature = extractor.extract(viccEntry);

        // TODO Should be multiple annotations
        assertEquals(1, codonsPerFeature.size());
        assertEquals(1, codonsPerFeature.values().size());
        CodonAnnotation annotation = codonsPerFeature.get(viccEntry.features().get(0)).get(0);

        assertEquals("12", annotation.chromosome());
        assertEquals(25380168, annotation.start());
        assertEquals(25378707, annotation.end());
        assertEquals("KRAS", annotation.gene());
        assertEquals(MutationTypeFilter.MISSENSE_ANY, annotation.mutationType());
    }

    @Test
    public void canExtractAnnotationForEntryWithTwoCodonFeatures() {
        CodonExtractor extractor = createWithDriverGenes(createDriverGenes("TP53", "PIK3CA", "KRAS"));

        Feature feature1 = ViccTestFactory.testFeatureWithGeneAndName("PIK3CA", "E545X");
        Feature feature2 = ViccTestFactory.testFeatureWithGeneAndName("KRAS", "G12X");
        ViccEntry entry = ViccTestFactory.testEntryWithFeatures(Lists.newArrayList(feature1, feature2));

        assertEquals(2, extractor.extract(entry).size());
        assertEquals(1, extractor.extract(entry).get(feature1).size());
        assertEquals(1, extractor.extract(entry).get(feature2).size());
    }

    @Test
    public void canExtractCodons() {
        assertEquals(Integer.valueOf("600"), CodonExtractor.extractCodonNumber("BRAF (V600)"));
        assertEquals(Integer.valueOf("742"), CodonExtractor.extractCodonNumber("W742"));
        assertEquals(Integer.valueOf("179"), CodonExtractor.extractCodonNumber("Q179X"));
        assertEquals(Integer.valueOf("61"), CodonExtractor.extractCodonNumber("KRAS Q61X"));

        assertNull(CodonExtractor.extractCodonNumber("Not a codon number"));
    }

    @NotNull
    private static CodonExtractor createWithDriverGenes(@NotNull List<DriverGene> driverGenes) {
        return new CodonExtractor(HG19_GENE_CHECKER, driverGenes, HG19_GENE_MAP);
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