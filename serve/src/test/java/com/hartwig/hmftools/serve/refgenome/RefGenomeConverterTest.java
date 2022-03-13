package com.hartwig.hmftools.serve.refgenome;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.google.common.collect.Sets;
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
import com.hartwig.hmftools.serve.actionability.range.RangeType;
import com.hartwig.hmftools.serve.extraction.codon.CodonAnnotation;
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
    public void canCurateKnownCodons() {
        Set<KnownCodon> knownCodonsSet = Sets.newHashSet();
        KnownCodon codon1 = ImmutableKnownCodon.builder()
                .from(ServeTestFactory.createTestKnownCodon())
                .annotation(ImmutableCodonAnnotation.builder()
                        .from(ServeTestFactory.createTestCodonAnnotation())
                        .gene("BRAF")
                        .chromosome(TEST_CHROMOSOME)
                        .start(10)
                        .end(20)
                        .rank(600)
                        .transcript("transcript")
                        .build())
                .build();
        KnownCodon codon2 = ImmutableKnownCodon.builder()
                .from(ServeTestFactory.createTestKnownCodon())
                .annotation(ImmutableCodonAnnotation.builder()
                        .from(ServeTestFactory.createTestCodonAnnotation())
                        .gene("BRAF")
                        .chromosome(TEST_CHROMOSOME)
                        .start(30)
                        .end(40)
                        .rank(601)
                        .transcript("transcript")
                        .build())
                .build();
        knownCodonsSet.add(codon1);
        knownCodonsSet.add(codon2);

        Set<KnownCodon> curatedCodon37 = RefGenomeConverter.curateKnownCodons(knownCodonsSet, RefGenomeVersion.V37);
        KnownCodon knownCodon1 = findByRankCodon(curatedCodon37, 601);
        assertEquals("BRAF", knownCodon1.annotation().gene());
        assertEquals(30, knownCodon1.annotation().start());
        assertEquals(40, knownCodon1.annotation().end());
        assertEquals("transcript", knownCodon1.annotation().transcript());

        KnownCodon knownCodon2 = findByRankCodon(curatedCodon37, 600);
        assertEquals("BRAF", knownCodon2.annotation().gene());
        assertEquals(140753335, knownCodon2.annotation().start());
        assertEquals(140753337, knownCodon2.annotation().end());
        assertEquals("ENST00000288602", knownCodon2.annotation().transcript());

        Set<KnownCodon> curatedCodon38 = RefGenomeConverter.curateKnownCodons(knownCodonsSet, RefGenomeVersion.V38);
        KnownCodon knownCodon3 = findByRankCodon(curatedCodon38, 601);
        assertEquals("BRAF", knownCodon3.annotation().gene());
        assertEquals(30, knownCodon3.annotation().start());
        assertEquals(40, knownCodon3.annotation().end());
        assertEquals("transcript", knownCodon3.annotation().transcript());

        KnownCodon knownCodon4 = findByRankCodon(curatedCodon38, 600);
        assertEquals("BRAF", knownCodon4.annotation().gene());
        assertEquals(10, knownCodon4.annotation().start());
        assertEquals(20, knownCodon4.annotation().end());
        assertEquals("transcript", knownCodon4.annotation().transcript());
    }

    @NotNull
    private static KnownCodon findByRankCodon(@NotNull Iterable<KnownCodon> knownCodons, int rank) {
        for (KnownCodon codon : knownCodons) {
            if (codon.annotation().rank() == rank) {
                return codon;
            }
        }

        throw new IllegalStateException("Could not find known codon with rank " + rank);
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
    public void canCurateActionableRanges() {
        Set<ActionableRange> actionableRangeSet = Sets.newHashSet();
        ActionableRange range1 = ImmutableActionableRange.builder()
                .from(ServeTestFactory.createTestActionableRange())
                .gene("BRAF")
                .rank(600)
                .rangeType(RangeType.CODON)
                .start(10)
                .end(20)
                .transcript("transcript")
                .build();
        ActionableRange range2 = ImmutableActionableRange.builder()
                .from(ServeTestFactory.createTestActionableRange())
                .gene("BRAF")
                .rank(601)
                .rangeType(RangeType.CODON)
                .start(30)
                .end(40)
                .transcript("transcript")
                .build();
        actionableRangeSet.add(range1);
        actionableRangeSet.add(range2);

        Set<ActionableRange> curatedRanges37 = RefGenomeConverter.curateActionableRange(actionableRangeSet, RefGenomeVersion.V37);
        ActionableRange actionableRange1 = findByRankRange(curatedRanges37, 601);
        assertEquals("BRAF", actionableRange1.gene());
        assertEquals(30, actionableRange1.start());
        assertEquals(40, actionableRange1.end());
        assertEquals("transcript", actionableRange1.transcript());

        ActionableRange actionableRange2 = findByRankRange(curatedRanges37, 600);
        assertEquals("BRAF", actionableRange2.gene());
        assertEquals(140753335, actionableRange2.start());
        assertEquals(140753337, actionableRange2.end());
        assertEquals("ENST00000288602", actionableRange2.transcript());

        Set<ActionableRange> curatedRanges38 = RefGenomeConverter.curateActionableRange(actionableRangeSet, RefGenomeVersion.V38);
        ActionableRange actionableRange3 = findByRankRange(curatedRanges38, 601);
        assertEquals("BRAF", actionableRange3.gene());
        assertEquals(30, actionableRange3.start());
        assertEquals(40, actionableRange3.end());
        assertEquals("transcript", actionableRange3.transcript());

        ActionableRange actionableRange4 = findByRankRange(curatedRanges38, 600);
        assertEquals("BRAF", actionableRange4.gene());
        assertEquals(10, actionableRange4.start());
        assertEquals(20, actionableRange4.end());
        assertEquals("transcript", actionableRange4.transcript());
    }

    @NotNull
    private static ActionableRange findByRankRange(@NotNull Iterable<ActionableRange> ranges, int rank) {
        for (ActionableRange range : ranges) {
            if (range.rank() == rank) {
                return range;
            }
        }

        throw new IllegalStateException("Could not find actionable range with rank " + rank);
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
                algo);
    }

    private static class DummyLiftOver implements LiftOverAlgo {

        @Nullable
        @Override
        public LiftOverResult liftOver(@NotNull final String chromosome, final int position) {
            return ImmutableLiftOverResult.builder().chromosome(chromosome).position(position).build();
        }
    }

    private static class NullLiftOver implements LiftOverAlgo {

        @Nullable
        @Override
        public LiftOverResult liftOver(@NotNull final String chromosome, final int position) {
            return null;
        }
    }
}