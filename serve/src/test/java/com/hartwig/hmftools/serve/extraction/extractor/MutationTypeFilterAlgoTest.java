package com.hartwig.hmftools.serve.extraction.extractor;

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

public class MutationTypeFilterAlgoTest {

    @Test
    public void canDetermineMutationFilter() {
        String tsg = "tsg";
        String onco = "onco";

        MutationTypeFilterAlgo algo = new MutationTypeFilterAlgo(createDriverGenes(tsg, onco));

        assertEquals(MutationTypeFilter.NONSENSE_OR_FRAMESHIFT, algo.determine(onco, "EXON 9 FRAMESHIFT"));
        assertEquals(MutationTypeFilter.SPLICE, algo.determine(onco, "Exon 12 splice site insertion"));
        assertEquals(MutationTypeFilter.SPLICE, algo.determine(onco, "EXON 14 SKIPPING MUTATION"));
        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_DELETION, algo.determine(onco, "EGFR exon 19 deletions"));
        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_DELETION, algo.determine(onco, "Null (Partial deletion of Exons 2 & 3)"));
        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_INSERTION, algo.determine(onco, "Exon 20 insertions"));
        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_ANY, algo.determine(onco, "Exon 20 insertions/deletions"));
        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_ANY, algo.determine(onco, "Exon 19 deletion/insertion"));
        assertEquals(MutationTypeFilter.MISSENSE_ANY, algo.determine(onco, "mut"));
        assertEquals(MutationTypeFilter.MISSENSE_INFRAME_INSERTION, algo.determine(onco, "insertion"));
        assertEquals(MutationTypeFilter.ANY, algo.determine(tsg, "mut"));
        assertEquals(MutationTypeFilter.NONSENSE_OR_FRAMESHIFT, algo.determine(tsg, "frameshift"));

        assertEquals(MutationTypeFilter.ANY, algo.determine("NOT-A-GENE", "abcd"));
    }

    @NotNull
    private static List<DriverGene> createDriverGenes(@NotNull String geneTsg, @NotNull String geneOnco) {
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
        DriverGene driverGeneOnco = driverGeneBuilder.gene(geneOnco).likelihoodType(ONCO).build();

        return Lists.newArrayList(driverGeneTsg, driverGeneOnco);
    }
}