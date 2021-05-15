package com.hartwig.hmftools.serve.refgenome;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.GeneNameMapping;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.serve.ServeTestFactory;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.fusion.ImmutableActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.gene.ImmutableActionableGene;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.hotspot.ImmutableActionableHotspot;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.actionability.range.ImmutableActionableRange;
import com.hartwig.hmftools.serve.extraction.codon.ImmutableCodonAnnotation;
import com.hartwig.hmftools.serve.extraction.codon.ImmutableKnownCodon;
import com.hartwig.hmftools.serve.extraction.codon.KnownCodon;
import com.hartwig.hmftools.serve.extraction.copynumber.ImmutableKnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.exon.ImmutableExonAnnotation;
import com.hartwig.hmftools.serve.extraction.exon.ImmutableKnownExon;
import com.hartwig.hmftools.serve.extraction.exon.KnownExon;
import com.hartwig.hmftools.serve.extraction.fusion.ImmutableKnownFusionPair;
import com.hartwig.hmftools.serve.extraction.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.extraction.hotspot.ImmutableKnownHotspot;
import com.hartwig.hmftools.serve.extraction.hotspot.KnownHotspot;
import com.hartwig.hmftools.serve.refgenome.liftover.ImmutableLiftOverResult;
import com.hartwig.hmftools.serve.refgenome.liftover.LiftOverAlgo;
import com.hartwig.hmftools.serve.refgenome.liftover.LiftOverResult;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class RefGenomeConverterTest {

    private static final RefGenomeConverter DUMMY_CONVERTER = build37To38DummyConverter();
    private static final RefGenomeConverter NULL_CONVERTER = build37To38NullConverter();

    private static final String TEST_GENE = "BRAF";
    private static final String TEST_CHROMOSOME = "chr1";

    @Test
    public void canMapNormalGenes() {
        Map<String, String> v37 = Maps.newHashMap();
        Map<String, String> v38 = Maps.newHashMap();

        v37.put("X", "Y");
        v38.put("Y", "X");

        RefGenomeConverter converter37To38 = new RefGenomeConverter(RefGenomeVersion.V37,
                RefGenomeVersion.V38,
                RefGenomeResourceTestFactory.loadTestRefSequence38(),
                new DummyLiftOver(),
                new GeneNameMapping(v37, v38));

        KnownCopyNumber copyNumber37 =
                ImmutableKnownCopyNumber.builder().from(ServeTestFactory.createTestKnownCopyNumber()).gene("X").build();

        Set<KnownCopyNumber> convertedCopyNumber38 = converter37To38.convertKnownCopyNumbers(Sets.newHashSet(copyNumber37));
        assertEquals("Y", convertedCopyNumber38.iterator().next().gene());

        RefGenomeConverter converter38To37 = new RefGenomeConverter(RefGenomeVersion.V38,
                RefGenomeVersion.V37,
                RefGenomeResourceTestFactory.loadTestRefSequence37(),
                new DummyLiftOver(),
                new GeneNameMapping(v37, v38));

        KnownCopyNumber copyNumber38 =
                ImmutableKnownCopyNumber.builder().from(ServeTestFactory.createTestKnownCopyNumber()).gene("Y").build();

        Set<KnownCopyNumber> convertedCopyNumber37 = converter38To37.convertKnownCopyNumbers(Sets.newHashSet(copyNumber38));
        assertEquals("X", convertedCopyNumber37.iterator().next().gene());
    }

    @Test
    public void canSkipInvalidGeneMapping() {
        RefGenomeConverter converter = build37To38NullConverter();
        KnownCopyNumber copyNumber38 =
                ImmutableKnownCopyNumber.builder().from(ServeTestFactory.createTestKnownCopyNumber()).gene("RYBP").build();

        assertNotNull(converter.convertKnownCopyNumbers(Sets.newHashSet(copyNumber38)));
    }

    @Test
    public void canConvertKnownHotspots() {
        KnownHotspot hotspot = ImmutableKnownHotspot.builder()
                .from(ServeTestFactory.createTestKnownHotspot())
                .gene(TEST_GENE)
                .chromosome(TEST_CHROMOSOME)
                .position(1)
                .ref("G")
                .alt("T")
                .addSources(Knowledgebase.HARTWIG_CURATED)
                .build();

        Set<KnownHotspot> convertedHotspots = DUMMY_CONVERTER.convertKnownHotspots(Sets.newHashSet(hotspot));
        assertEquals(hotspot, convertedHotspots.iterator().next());

        assertTrue(NULL_CONVERTER.convertKnownHotspots(Sets.newHashSet(hotspot)).isEmpty());
    }

    @Test
    public void canConvertKnownCodons() {
        KnownCodon codon = ImmutableKnownCodon.builder()
                .from(ServeTestFactory.createTestKnownCodon())
                .annotation(ImmutableCodonAnnotation.builder()
                        .from(ServeTestFactory.createTestCodonAnnotation())
                        .gene(TEST_GENE)
                        .chromosome(TEST_CHROMOSOME)
                        .start(1)
                        .end(3)
                        .build())
                .addSources(Knowledgebase.HARTWIG_CURATED)
                .build();

        Set<KnownCodon> convertedCodons = DUMMY_CONVERTER.convertKnownCodons(Sets.newHashSet(codon));
        assertEquals(codon, convertedCodons.iterator().next());

        assertTrue(NULL_CONVERTER.convertKnownCodons(Sets.newHashSet(codon)).isEmpty());
    }

    @Test
    public void canConvertKnownExons() {
        KnownExon exon = ImmutableKnownExon.builder()
                .from(ServeTestFactory.createTestKnownExon())
                .annotation(ImmutableExonAnnotation.builder()
                        .from(ServeTestFactory.createTestExonAnnotation())
                        .gene(TEST_GENE)
                        .chromosome(TEST_CHROMOSOME)
                        .start(1)
                        .end(7)
                        .build())
                .addSources(Knowledgebase.HARTWIG_CURATED)
                .build();

        Set<KnownExon> convertedExons = DUMMY_CONVERTER.convertKnownExons(Sets.newHashSet(exon));
        assertEquals(exon, convertedExons.iterator().next());

        assertTrue(NULL_CONVERTER.convertKnownExons(Sets.newHashSet(exon)).isEmpty());
    }

    @Test
    public void canConvertKnownCopyNumbers() {
        KnownCopyNumber copyNumber = ImmutableKnownCopyNumber.builder()
                .from(ServeTestFactory.createTestKnownCopyNumber())
                .gene(TEST_GENE)
                .addSources(Knowledgebase.HARTWIG_CURATED)
                .build();

        Set<KnownCopyNumber> convertedCopyNumbers = DUMMY_CONVERTER.convertKnownCopyNumbers(Sets.newHashSet(copyNumber));
        assertEquals(copyNumber, convertedCopyNumbers.iterator().next());
    }

    @Test
    public void canConvertKnownFusionPairs() {
        KnownFusionPair fusionPair = ImmutableKnownFusionPair.builder()
                .from(ServeTestFactory.createTestKnownFusionPair())
                .geneUp(TEST_GENE)
                .geneDown(TEST_GENE)
                .addSources(Knowledgebase.HARTWIG_CURATED)
                .build();

        Set<KnownFusionPair> convertedFusionPairs = DUMMY_CONVERTER.convertKnownFusionPairs(Sets.newHashSet(fusionPair));
        assertEquals(fusionPair, convertedFusionPairs.iterator().next());
    }

    @Test
    public void canConvertActionableHotspots() {
        ActionableHotspot actionableHotspot = ImmutableActionableHotspot.builder()
                .from(ServeTestFactory.createTestActionableHotspotForSource(Knowledgebase.HARTWIG_CURATED))
                .chromosome(TEST_CHROMOSOME)
                .position(1)
                .ref("G")
                .alt("C")
                .build();

        Set<ActionableHotspot> convertedActionableHotspots = DUMMY_CONVERTER.convertActionableHotspots(Sets.newHashSet(actionableHotspot));
        assertEquals(actionableHotspot, convertedActionableHotspots.iterator().next());
    }

    @Test
    public void canConvertActionableRanges() {
        ActionableRange actionableRange = ImmutableActionableRange.builder()
                .from(ServeTestFactory.createTestActionableRangeForSource(Knowledgebase.HARTWIG_CURATED))
                .gene(TEST_GENE)
                .chromosome(TEST_CHROMOSOME)
                .start(3)
                .end(4)
                .build();

        Set<ActionableRange> convertedActionableRanges = DUMMY_CONVERTER.convertActionableRanges(Sets.newHashSet(actionableRange));
        assertEquals(actionableRange, convertedActionableRanges.iterator().next());
    }

    @Test
    public void canConvertActionableGenes() {
        ActionableGene actionableGene = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGeneForSource(Knowledgebase.HARTWIG_CURATED))
                .gene(TEST_GENE)
                .build();

        Set<ActionableGene> convertedActionableGenes = DUMMY_CONVERTER.convertActionableGenes(Sets.newHashSet(actionableGene));
        assertEquals(actionableGene, convertedActionableGenes.iterator().next());
    }

    @Test
    public void canConvertActionableFusions() {
        ActionableFusion actionableFusion = ImmutableActionableFusion.builder()
                .from(ServeTestFactory.createTestActionableFusionForSource(Knowledgebase.HARTWIG_CURATED))
                .geneUp(TEST_GENE)
                .geneDown(TEST_GENE)
                .build();

        Set<ActionableFusion> convertedActionableFusions = DUMMY_CONVERTER.convertActionableFusions(Sets.newHashSet(actionableFusion));
        assertEquals(actionableFusion, convertedActionableFusions.iterator().next());
    }

    @NotNull
    private static RefGenomeConverter build37To38DummyConverter() {
        return build37To38ConverterWithLiftOverAlgo(new DummyLiftOver());
    }

    @NotNull
    private static RefGenomeConverter build37To38NullConverter() {
        return build37To38ConverterWithLiftOverAlgo(new NullLiftOver());
    }

    @NotNull
    private static RefGenomeConverter build37To38ConverterWithLiftOverAlgo(@NotNull LiftOverAlgo algo) {
        return new RefGenomeConverter(RefGenomeVersion.V37,
                RefGenomeVersion.V38,
                RefGenomeResourceTestFactory.loadTestRefSequence38(),
                algo,
                GeneNameMapping.loadFromEmbeddedResource());
    }

    private static class DummyLiftOver implements LiftOverAlgo {

        @Nullable
        @Override
        public LiftOverResult liftOver(@NotNull final String chromosome, final long position) {
            return ImmutableLiftOverResult.builder().chromosome(chromosome).position(position).build();
        }
    }

    private static class NullLiftOver implements LiftOverAlgo {

        @Nullable
        @Override
        public LiftOverResult liftOver(@NotNull final String chromosome, final long position) {
            return null;
        }
    }
}