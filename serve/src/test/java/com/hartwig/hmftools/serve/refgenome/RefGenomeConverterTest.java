package com.hartwig.hmftools.serve.refgenome;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.GeneNameMapping;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.serve.ServeTestFactory;
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

        // Codons that are not 3 bases long can't get converted
        KnownCodon invalidCodon = ImmutableKnownCodon.builder()
                .from(ServeTestFactory.createTestKnownCodon())
                .annotation(ImmutableCodonAnnotation.builder()
                        .from(ServeTestFactory.createTestCodonAnnotation())
                        .gene(TEST_GENE)
                        .chromosome(TEST_CHROMOSOME)
                        .start(1)
                        .end(2)
                        .build())
                .addSources(Knowledgebase.HARTWIG_CURATED)
                .build();

        assertTrue(DUMMY_CONVERTER.convertKnownCodons(Sets.newHashSet(invalidCodon)).isEmpty());
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