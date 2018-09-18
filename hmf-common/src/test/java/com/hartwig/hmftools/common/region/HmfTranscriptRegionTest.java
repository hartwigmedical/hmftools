package com.hartwig.hmftools.common.region;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class HmfTranscriptRegionTest {

    //MIVO:     both forward and reverse cases use the following scenario (* = UTR):
    //exons:        |--------------------|   |------------------------------------------------|
    //positions:    --1*---2---3---4---5---x---7---8---9---10---11---12---13---14---15*---16*--
    //codons:            |___________|_______________|_____________|______________|

    @Test
    public void canRetrieveExonByIndexOnForwardGene() {
        HmfTranscriptRegion transcript = create(Strand.FORWARD);

        assertNull(transcript.exonByIndex(0));

        HmfExonRegion exon1 = transcript.exonByIndex(1);
        assertNotNull(exon1);
        assertEquals("1", exon1.exonID());

        HmfExonRegion exon2 = transcript.exonByIndex(2);
        assertNotNull(exon2);
        assertEquals("2", exon2.exonID());

        assertNull(transcript.exonByIndex(3));
    }

    @Test
    public void canRetrieveExonByIndexOnReverseGene() {
        HmfTranscriptRegion transcript = create(Strand.REVERSE);

        assertNull(transcript.exonByIndex(0));

        HmfExonRegion exon1 = transcript.exonByIndex(1);
        assertNotNull(exon1);
        assertEquals("2", exon1.exonID());

        HmfExonRegion exon2 = transcript.exonByIndex(2);
        assertNotNull(exon2);
        assertEquals("1", exon2.exonID());

        assertNull(transcript.exonByIndex(3));
    }

    @Test
    public void canRetrieveCodonByIndexOnForwardGene() {
        HmfTranscriptRegion transcript = create(Strand.FORWARD);

        assertNull(transcript.codonByIndex(-1));
        assertNull(transcript.codonByIndex(0));

        GenomeRegion codon1 = assertedCodonGet(transcript, 1);
        assertNotNull(codon1);
        assertEquals(2, codon1.start());
        assertEquals(4, codon1.end());

        GenomeRegion codon2Part1 = assertedCodonGet(transcript, 2, 0);
        assertNotNull(codon2Part1);
        assertEquals(5, codon2Part1.start());
        assertEquals(5, codon2Part1.end());

        GenomeRegion codon2Part2 = assertedCodonGet(transcript, 2, 1);
        assertNotNull(codon2Part2);
        assertEquals(7, codon2Part2.start());
        assertEquals(8, codon2Part2.end());

        GenomeRegion codon3 = assertedCodonGet(transcript, 3);
        assertNotNull(codon3);
        assertEquals(9, codon3.start());
        assertEquals(11, codon3.end());

        GenomeRegion codon4 = assertedCodonGet(transcript, 4);
        assertNotNull(codon4);
        assertEquals(12, codon4.start());
        assertEquals(14, codon4.end());

        assertNull(transcript.codonByIndex(5));
    }

    @Test
    public void canRetrieveCodonByIndexOnReverseGene() {
        HmfTranscriptRegion transcript = create(Strand.REVERSE);

        assertNull(transcript.codonByIndex(-1));
        assertNull(transcript.codonByIndex(0));

        GenomeRegion codon1 = assertedCodonGet(transcript, 1);
        assertNotNull(codon1);
        assertEquals(12, codon1.start());
        assertEquals(14, codon1.end());

        GenomeRegion codon2 = assertedCodonGet(transcript, 2);
        assertNotNull(codon2);
        assertEquals(9, codon2.start());
        assertEquals(11, codon2.end());

        GenomeRegion codon3Part1 = assertedCodonGet(transcript, 3, 0);
        assertNotNull(codon3Part1);
        assertEquals(5, codon3Part1.start());
        assertEquals(5, codon3Part1.end());

        GenomeRegion codon3Part2 = assertedCodonGet(transcript, 3, 1);
        assertNotNull(codon3Part2);
        assertEquals(7, codon3Part2.start());
        assertEquals(8, codon3Part2.end());

        GenomeRegion codon4 = assertedCodonGet(transcript, 4);
        assertNotNull(codon4);
        assertEquals(2, codon4.start());
        assertEquals(4, codon4.end());

        assertNull(transcript.codonByIndex(5));
    }

    @Test
    public void canRetrieveCodonRangeByIndexOnForwardGene() {
        HmfTranscriptRegion transcript = create(Strand.FORWARD);

        assertNull(transcript.codonRangeByIndex(-1, 1));
        assertNull(transcript.codonRangeByIndex(0, 0));

        List<GenomeRegion> codonsOneAndTwo = transcript.codonRangeByIndex(1, 2);
        assertNotNull(codonsOneAndTwo);
        GenomeRegion codonsOneAndTwoPart1 = codonsOneAndTwo.get(0);
        assertEquals(2, codonsOneAndTwoPart1.start());
        assertEquals(5, codonsOneAndTwoPart1.end());

        GenomeRegion codonsOneAndTwoPart2 = codonsOneAndTwo.get(1);
        assertEquals(7, codonsOneAndTwoPart2.start());
        assertEquals(8, codonsOneAndTwoPart2.end());

        List<GenomeRegion> allCodons = transcript.codonRangeByIndex(1, 4);
        assertNotNull(allCodons);
        GenomeRegion allCodonsPart1 = allCodons.get(0);
        assertEquals(2, allCodonsPart1.start());
        assertEquals(5, allCodonsPart1.end());

        GenomeRegion allCodonsPart2 = allCodons.get(1);
        assertEquals(7, allCodonsPart2.start());
        assertEquals(14, allCodonsPart2.end());

        assertNull(transcript.codonRangeByIndex(0, 5));
    }

    @Test
    public void canRetrieveCodonRangeByIndexOnReverseGene() {
        HmfTranscriptRegion transcript = create(Strand.REVERSE);

        assertNull(transcript.codonRangeByIndex(-1, 1));
        assertNull(transcript.codonRangeByIndex(0, 0));

        List<GenomeRegion> codonsOneAndTwo = transcript.codonRangeByIndex(1, 2);
        assertNotNull(codonsOneAndTwo);
        assertEquals(1, codonsOneAndTwo.size());
        GenomeRegion codonsOneAndTwoPart1 = codonsOneAndTwo.get(0);
        assertEquals(9, codonsOneAndTwoPart1.start());
        assertEquals(14, codonsOneAndTwoPart1.end());

        List<GenomeRegion> allCodons = transcript.codonRangeByIndex(1, 4);
        assertNotNull(allCodons);
        GenomeRegion allCodonsPart1 = allCodons.get(0);
        assertEquals(2, allCodonsPart1.start());
        assertEquals(5, allCodonsPart1.end());

        GenomeRegion allCodonsPart2 = allCodons.get(1);
        assertEquals(7, allCodonsPart2.start());
        assertEquals(14, allCodonsPart2.end());

        assertNull(transcript.codonRangeByIndex(0, 5));
    }

    @Test
    public void canRetrieveCodingRangeByGenomicCoordinatesOnForwardGene() {
        runCodingRangeByGenomicCoordinatesTest(create(Strand.FORWARD));
    }

    @Test
    public void canRetrieveCodingRangeByGenomicCoordinatesOnReverseGene() {
        runCodingRangeByGenomicCoordinatesTest(create(Strand.REVERSE));
    }

    private static void runCodingRangeByGenomicCoordinatesTest(@NotNull HmfTranscriptRegion transcript) {
        assertNull(transcript.codingRangeByGenomicCoordinates(-1, 1));
        assertNull(transcript.codingRangeByGenomicCoordinates(0, 0));

        List<GenomeRegion> codingRangePartial = transcript.codingRangeByGenomicCoordinates(1, 12);
        assertNotNull(codingRangePartial);
        GenomeRegion codingRangePartial1 = codingRangePartial.get(0);
        assertEquals(2, codingRangePartial1.start());
        assertEquals(5, codingRangePartial1.end());

        GenomeRegion codingRangePartial2 = codingRangePartial.get(1);
        assertEquals(7, codingRangePartial2.start());
        assertEquals(12, codingRangePartial2.end());

        assertNull(transcript.codingRangeByGenomicCoordinates(0, 5));
    }

    @Test
    public void worksForRealBRAFCodon600() {
        HmfTranscriptRegion braf = HmfGenePanelSupplier.allGenesMap().get("BRAF");

        List<GenomeRegion> codon600 = braf.codonByIndex(600);
        assertNotNull(codon600);
        assertEquals(1, codon600.size());

        assertEquals(140453135, codon600.get(0).start());
        assertEquals(140453137, codon600.get(0).end());
    }

    @NotNull
    private static GenomeRegion assertedCodonGet(@NotNull HmfTranscriptRegion transcript, int codonIndex) {
        return assertedCodonGet(transcript, codonIndex, 0);
    }

    @NotNull
    private static GenomeRegion assertedCodonGet(@NotNull HmfTranscriptRegion transcript, int codonIndex, int regionIndex) {
        List<GenomeRegion> regions = transcript.codonByIndex(codonIndex);
        assertNotNull(regions);
        return regions.get(regionIndex);
    }

    @NotNull
    private static HmfTranscriptRegion create(@NotNull Strand strand) {
        List<HmfExonRegion> exome = Lists.newArrayList(ImmutableHmfExonRegion.builder().chromosome("1").exonID("1").start(1).end(5).build(),
                ImmutableHmfExonRegion.builder().chromosome("1").exonID("2").start(7).end(16).build());

        return ImmutableHmfTranscriptRegion.builder()
                .chromosome("1")
                .chromosomeBand("band")
                .gene("gene")
                .geneID("geneID")
                .transcriptID("transcriptID")
                .transcriptVersion(1)
                .geneStart(1)
                .geneEnd(16)
                .exome(exome)
                .start(1)
                .end(16)
                .codingStart(2)
                .codingEnd(14)
                .strand(strand)
                .build();

    }
}