package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.ONCO;
import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.TSG;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;
import com.hartwig.hmftools.serve.actionability.range.MutationTypeFilter;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class MutationTypeFilterExtractionTest {

    @Test
    public void canExtractMutationFilter() {
        List<DriverGene> driverGenes = createDriverGenes("TP53", "EGFR", "ERBB2");
        String gene = "ERBB2";

        assertEquals(MutationTypeFilter.NONSENSE_OR_FRAMESHIFT,
                MutationTypeFilterExtraction.extract("EXON 9 FRAMESHIFT", driverGenes, gene));

        assertEquals(MutationTypeFilter.SPLICE, MutationTypeFilterExtraction.extract("Exon 12 splice site insertion", driverGenes, gene));

        assertEquals(MutationTypeFilter.SPLICE, MutationTypeFilterExtraction.extract("EXON 14 SKIPPING MUTATION", driverGenes, gene));

        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_DELETION,
                MutationTypeFilterExtraction.extract("EGFR exon 19 deletions", driverGenes, gene));

        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_DELETION,
                MutationTypeFilterExtraction.extract("Null (Partial deletion of Exons 2 & 3)", driverGenes, gene));

        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_INSERTION,
                MutationTypeFilterExtraction.extract("Exon 20 insertions", driverGenes, gene));

        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_ANY,
                MutationTypeFilterExtraction.extract("Exon 20 insertions/deletions", driverGenes, gene));

        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_ANY,
                MutationTypeFilterExtraction.extract("Exon 19 deletion/insertion", driverGenes, gene));

        assertEquals(MutationTypeFilter.UNKNOWN, MutationTypeFilterExtraction.extract("abcd", driverGenes, "efgh"));

        assertEquals(MutationTypeFilter.MISSENSE_ANY, MutationTypeFilterExtraction.extract("mut", driverGenes, gene));

        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_INSERTION, MutationTypeFilterExtraction.extract("insertion", driverGenes, gene));

        assertEquals(MutationTypeFilter.ANY, MutationTypeFilterExtraction.extract("mut", driverGenes, "TP53"));

        assertEquals(MutationTypeFilter.NONSENSE_OR_FRAMESHIFT, MutationTypeFilterExtraction.extract("frameshift", driverGenes, "TP53"));
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