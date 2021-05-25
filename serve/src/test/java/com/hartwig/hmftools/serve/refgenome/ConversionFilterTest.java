package com.hartwig.hmftools.serve.refgenome;

import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.serve.ServeTestFactory;
import com.hartwig.hmftools.serve.actionability.fusion.ImmutableActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ImmutableActionableGene;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.ImmutableExtractionResult;
import com.hartwig.hmftools.serve.extraction.codon.ImmutableCodonAnnotation;
import com.hartwig.hmftools.serve.extraction.codon.ImmutableKnownCodon;
import com.hartwig.hmftools.serve.extraction.copynumber.ImmutableKnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.exon.ImmutableExonAnnotation;
import com.hartwig.hmftools.serve.extraction.exon.ImmutableKnownExon;
import com.hartwig.hmftools.serve.extraction.fusion.ImmutableKnownFusionPair;
import com.hartwig.hmftools.serve.extraction.hotspot.ImmutableKnownHotspot;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ConversionFilterTest {

    @Test
    public void canFilterGenes() {
        ConversionFilter filter = new ConversionFilter();

        ExtractionResult resultToFilter = createExtractionResultForGene(ConversionFilterFactory.GENES_TO_FILTER.iterator().next());
        ExtractionResult filtered = filter.filter(resultToFilter);
        assertTrue(filtered.knownHotspots().isEmpty());
        assertTrue(filtered.knownCodons().isEmpty());
        assertTrue(filtered.knownExons().isEmpty());
        assertTrue(filtered.knownCopyNumbers().isEmpty());
        assertTrue(filtered.knownFusionPairs().isEmpty());
        assertTrue(filtered.actionableGenes().isEmpty());
        assertTrue(filtered.actionableFusions().isEmpty());

        filter.reportUnusedFilterEntries();
    }

    @NotNull
    private static ExtractionResult createExtractionResultForGene(@NotNull String gene) {
        return ImmutableExtractionResult.builder()
                .refGenomeVersion(RefGenomeVersion.V38)
                .addKnownHotspots(ImmutableKnownHotspot.builder().from(ServeTestFactory.createTestKnownHotspot()).gene(gene).build())
                .addKnownCodons(ImmutableKnownCodon.builder()
                        .from(ServeTestFactory.createTestKnownCodon())
                        .annotation(ImmutableCodonAnnotation.builder()
                                .from(ServeTestFactory.createTestCodonAnnotation())
                                .gene(gene)
                                .build())
                        .build())
                .addKnownExons(ImmutableKnownExon.builder()
                        .from(ServeTestFactory.createTestKnownExon())
                        .annotation(ImmutableExonAnnotation.builder().from(ServeTestFactory.createTestExonAnnotation()).gene(gene).build())
                        .build())
                .addKnownCopyNumbers(ImmutableKnownCopyNumber.builder()
                        .from(ServeTestFactory.createTestKnownCopyNumber())
                        .gene(gene)
                        .build())
                .addKnownFusionPairs(ImmutableKnownFusionPair.builder()
                        .from(ServeTestFactory.createTestKnownFusionPair())
                        .geneUp(gene)
                        .build())
                .addKnownFusionPairs(ImmutableKnownFusionPair.builder()
                        .from(ServeTestFactory.createTestKnownFusionPair())
                        .geneDown(gene)
                        .build())
                .addActionableGenes(ImmutableActionableGene.builder().from(ServeTestFactory.createTestActionableGene()).gene(gene).build())
                .addActionableFusions(ImmutableActionableFusion.builder()
                        .from(ServeTestFactory.createTestActionableFusion())
                        .geneUp(gene)
                        .build())
                .addActionableFusions(ImmutableActionableFusion.builder()
                        .from(ServeTestFactory.createTestActionableFusion())
                        .geneDown(gene)
                        .build())
                .build();
    }

}