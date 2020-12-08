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
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.actionability.range.MutationTypeFilter;
import com.hartwig.hmftools.serve.codon.CodonAnnotation;
import com.hartwig.hmftools.serve.extraction.GeneChecker;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CodonExtractorTest {

    private static final Map<String, HmfTranscriptRegion> HG19_GENE_MAP = HmfGenePanelSupplier.allGenesMap37();
    private static final GeneChecker HG19_GENE_CHECKER = new GeneChecker(HG19_GENE_MAP.keySet());

    @Test
    public void canExtractSimpleCodon() {
        CodonExtractor extractor = createWithDriverGenes(createDriverGenes("TP53", "EGFR", "KIT"));
        List<CodonAnnotation> codons = extractor.extract("TP53", null, EventType.CODON, "R249");

        assertEquals(1, codons.size());

        assertEquals("17", codons.get(0).chromosome());
        assertEquals(7577534, codons.get(0).start());
        assertEquals(7577536, codons.get(0).end());
        assertEquals("TP53", codons.get(0).gene());
        assertEquals(MutationTypeFilter.ANY, codons.get(0).mutationType());
    }

    @Test
    public void canExtractCodonOnMultipleExons() {
        CodonExtractor extractor = createWithDriverGenes(createDriverGenes("TP53", "KRAS", "KIT"));
        List<CodonAnnotation> codons = extractor.extract("KRAS", null, EventType.CODON, "R97");

        // TODO Should be multiple annotations
        assertEquals(1, codons.size());

        assertEquals("12", codons.get(0).chromosome());
        assertEquals(25380168, codons.get(0).start());
        assertEquals(25378707, codons.get(0).end());
        assertEquals("KRAS", codons.get(0).gene());
        assertEquals(MutationTypeFilter.MISSENSE_ANY, codons.get(0).mutationType());
    }

    @Test
    public void canExtractCodons() {
        assertEquals(600, (int) CodonExtractor.extractCodonNumber("BRAF (V600)"));
        assertEquals(742, (int) CodonExtractor.extractCodonNumber("W742"));
        assertEquals(179, (int) CodonExtractor.extractCodonNumber("Q179X"));
        assertEquals(61, (int) CodonExtractor.extractCodonNumber("KRAS Q61X"));

        assertNull(CodonExtractor.extractCodonNumber("Not a codon number"));
    }

    @NotNull
    private static CodonExtractor createWithDriverGenes(@NotNull List<DriverGene> driverGenes) {
        return new CodonExtractor(HG19_GENE_CHECKER, new MutationTypeFilterAlgo(driverGenes), HG19_GENE_MAP);
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