package com.hartwig.hmftools.common.region;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class HmfTranscriptRegionTest {

    //MIVO:     both forward and reverse cases use the following scenario (* = UTR):
    //exons:        |--------------------|   |------------------------------------------------|
    //positions:    --1*---2---3---4---5---x---7---8---9---10---11---12---13---14---15*---16*--
    //codons:       |____________|_______________|____________|______________|________________|

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
        assertEquals(1, codon1.start());
        assertEquals(3, codon1.end());

        GenomeRegion codon2Part1 = assertedCodonGet(transcript, 2, 0);
        assertNotNull(codon2Part1);
        assertEquals(4, codon2Part1.start());
        assertEquals(5, codon2Part1.end());

        GenomeRegion codon2Part2 = assertedCodonGet(transcript, 2, 1);
        assertNotNull(codon2Part2);
        assertEquals(7, codon2Part2.start());
        assertEquals(7, codon2Part2.end());

        GenomeRegion codon3 = assertedCodonGet(transcript, 3);
        assertNotNull(codon3);
        assertEquals(8, codon3.start());
        assertEquals(10, codon3.end());

        GenomeRegion codon4 = assertedCodonGet(transcript, 4);
        assertNotNull(codon4);
        assertEquals(11, codon4.start());
        assertEquals(13, codon4.end());

        GenomeRegion codon5 = assertedCodonGet(transcript, 5);
        assertNotNull(codon5);
        assertEquals(14, codon5.start());
        assertEquals(16, codon5.end());

        assertNull(transcript.codonByIndex(6));
    }

    @Test
    public void canRetrieveCodonByIndexOnReverseGene() {
        HmfTranscriptRegion transcript = create(Strand.REVERSE);

        assertNull(transcript.codonByIndex(-1));
        assertNull(transcript.codonByIndex(0));

        GenomeRegion codon1 = assertedCodonGet(transcript, 1);
        assertNotNull(codon1);
        assertEquals(14, codon1.start());
        assertEquals(16, codon1.end());

        GenomeRegion codon2 = assertedCodonGet(transcript, 2);
        assertNotNull(codon2);
        assertEquals(11, codon2.start());
        assertEquals(13, codon2.end());

        GenomeRegion codon3 = assertedCodonGet(transcript, 3);
        assertNotNull(codon3);
        assertEquals(8, codon3.start());
        assertEquals(10, codon3.end());

        GenomeRegion codon4Part1 = assertedCodonGet(transcript, 4, 0);
        assertNotNull(codon4Part1);
        assertEquals(4, codon4Part1.start());
        assertEquals(5, codon4Part1.end());

        GenomeRegion codon4Part2 = assertedCodonGet(transcript, 4, 1);
        assertNotNull(codon4Part2);
        assertEquals(7, codon4Part2.start());
        assertEquals(7, codon4Part2.end());

        GenomeRegion codon5 = assertedCodonGet(transcript, 5);
        assertNotNull(codon5);
        assertEquals(1, codon5.start());
        assertEquals(3, codon5.end());

        assertNull(transcript.codonByIndex(6));
    }

    @Test
    public void canRetrieveCodonRangeByIndexOnForwardGene() {
        HmfTranscriptRegion transcript = create(Strand.FORWARD);

        assertNull(transcript.codonRangeByIndex(-1, 1));
        assertNull(transcript.codonRangeByIndex(0, 0));

        List<GenomeRegion> codonsOneAndTwo = transcript.codonRangeByIndex(1, 2);
        assertNotNull(codonsOneAndTwo);
        GenomeRegion codonsOneAndTwoPart1 = codonsOneAndTwo.get(0);
        assertEquals(1, codonsOneAndTwoPart1.start());
        assertEquals(5, codonsOneAndTwoPart1.end());

        GenomeRegion codonsOneAndTwoPart2 = codonsOneAndTwo.get(1);
        assertEquals(7, codonsOneAndTwoPart2.start());
        assertEquals(7, codonsOneAndTwoPart2.end());

        List<GenomeRegion> allCodons = transcript.codonRangeByIndex(1, 5);
        assertNotNull(allCodons);
        GenomeRegion allCodonsPart1 = allCodons.get(0);
        assertEquals(1, allCodonsPart1.start());
        assertEquals(5, allCodonsPart1.end());

        GenomeRegion allCodonsPart2 = allCodons.get(1);
        assertEquals(7, allCodonsPart2.start());
        assertEquals(16, allCodonsPart2.end());

        assertNull(transcript.codonRangeByIndex(0, 6));
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
        assertEquals(11, codonsOneAndTwoPart1.start());
        assertEquals(16, codonsOneAndTwoPart1.end());

        List<GenomeRegion> allCodons = transcript.codonRangeByIndex(1, 5);
        assertNotNull(allCodons);
        GenomeRegion allCodonsPart1 = allCodons.get(0);
        assertEquals(1, allCodonsPart1.start());
        assertEquals(5, allCodonsPart1.end());

        GenomeRegion allCodonsPart2 = allCodons.get(1);
        assertEquals(7, allCodonsPart2.start());
        assertEquals(16, allCodonsPart2.end());

        assertNull(transcript.codonRangeByIndex(0, 6));
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
                .codingStart(2)
                .start(1)
                .end(16)
                .codingEnd(14)
                .strand(strand)
                .build();

    }
}