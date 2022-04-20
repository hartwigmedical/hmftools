package com.hartwig.hmftools.common.genome.region;

import static com.hartwig.hmftools.common.genome.region.HmfTranscriptRegionUtils.codonByIndex;
import static com.hartwig.hmftools.common.genome.region.HmfTranscriptRegionUtils.codonRangeByRank;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CodonRegionTest
{

    // Both forward and reverse cases use the following scenario (* = UTR):
    // exons:        |--------------------|   |----------------------------------------------|   |------------------|
    // positions:    --1*---2---3---4---5---x---7---8---9---10---11---12---13---14---15---16---x---18---19*---20*--
    // codons:            |___________|_______________|_____________|______________|__________________|

    @Test
    public void canRetrieveExonByIndexOnForwardGene() {
        HmfTranscriptRegion transcript = create(Strand.FORWARD);

        assertNull(transcript.exonByIndex(0));

        HmfExonRegion exon1 = transcript.exonByIndex(1);
        assertNotNull(exon1);
        assertEquals(1, exon1.exonRank());

        HmfExonRegion exon2 = transcript.exonByIndex(2);
        assertNotNull(exon2);
        assertEquals(2, exon2.exonRank());

        HmfExonRegion exon3 = transcript.exonByIndex(3);
        assertNotNull(exon3);
        assertEquals(3, exon3.exonRank());

        assertNull(transcript.exonByIndex(4));
    }

    @Test
    public void canRetrieveExonByIndexOnReverseGene() {
        HmfTranscriptRegion transcript = create(Strand.REVERSE);

        assertNull(transcript.exonByIndex(0));

        HmfExonRegion exon1 = transcript.exonByIndex(1);
        assertNotNull(exon1);
        assertEquals(3, exon1.exonRank());

        HmfExonRegion exon2 = transcript.exonByIndex(2);
        assertNotNull(exon2);
        assertEquals(2, exon2.exonRank());

        HmfExonRegion exon3 = transcript.exonByIndex(3);
        assertNotNull(exon3);
        assertEquals(1, exon3.exonRank());

        assertNull(transcript.exonByIndex(4));
    }

    @Test
    public void canRetrieveCodonByIndexOnForwardGene() {
        HmfTranscriptRegion transcript = create(Strand.FORWARD);

        assertNull(codonByIndex(transcript, -1));
        assertNull(codonByIndex(transcript, 0));

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

        GenomeRegion codon5Part1 = assertedCodonGet(transcript, 5, 0);
        assertNotNull(codon5Part1);
        assertEquals(15, codon5Part1.start());
        assertEquals(16, codon5Part1.end());

        GenomeRegion codon5Part2 = assertedCodonGet(transcript, 5, 1);
        assertNotNull(codon5Part2);
        assertEquals(18, codon5Part2.start());
        assertEquals(18, codon5Part2.end());

        assertNull(codonByIndex(transcript, 6));
    }

    @Test
    public void canRetrieveCodonByIndexOnReverseGene() {
        HmfTranscriptRegion transcript = create(Strand.REVERSE);

        assertNull(codonByIndex(transcript, -1));
        assertNull(codonByIndex(transcript, 0));

        GenomeRegion codon1Part1 = assertedCodonGet(transcript, 1, 0);
        assertNotNull(codon1Part1);
        assertEquals(15, codon1Part1.start());
        assertEquals(16, codon1Part1.end());

        GenomeRegion codon1Part2 = assertedCodonGet(transcript, 1, 1);
        assertNotNull(codon1Part2);
        assertEquals(18, codon1Part2.start());
        assertEquals(18, codon1Part2.end());

        GenomeRegion codon2 = assertedCodonGet(transcript, 2);
        assertNotNull(codon2);
        assertEquals(12, codon2.start());
        assertEquals(14, codon2.end());

        GenomeRegion codon3 = assertedCodonGet(transcript, 3);
        assertNotNull(codon3);
        assertEquals(9, codon3.start());
        assertEquals(11, codon3.end());

        GenomeRegion codon4Part1 = assertedCodonGet(transcript, 4, 0);
        assertNotNull(codon4Part1);
        assertEquals(5, codon4Part1.start());
        assertEquals(5, codon4Part1.end());

        GenomeRegion codon4Part2 = assertedCodonGet(transcript, 4, 1);
        assertNotNull(codon4Part2);
        assertEquals(7, codon4Part2.start());
        assertEquals(8, codon4Part2.end());

        GenomeRegion codon5 = assertedCodonGet(transcript, 5);
        assertNotNull(codon5);
        assertEquals(2, codon5.start());
        assertEquals(4, codon5.end());

        assertNull(codonByIndex(transcript, 6));
    }

    @Test
    public void canRetrieveCodonRangeByIndexOnForwardGene() {
        HmfTranscriptRegion transcript = create(Strand.FORWARD);

        assertNull(codonRangeByRank(transcript, -1, 1));
        assertNull(codonRangeByRank(transcript, 0, 0));

        List<GenomeRegion> codonsOneAndTwo = codonRangeByRank(transcript, 1, 2);
        assertNotNull(codonsOneAndTwo);
        GenomeRegion codonsOneAndTwoPart1 = codonsOneAndTwo.get(0);
        assertEquals(2, codonsOneAndTwoPart1.start());
        assertEquals(5, codonsOneAndTwoPart1.end());

        GenomeRegion codonsOneAndTwoPart2 = codonsOneAndTwo.get(1);
        assertEquals(7, codonsOneAndTwoPart2.start());
        assertEquals(8, codonsOneAndTwoPart2.end());

        List<GenomeRegion> allCodons = codonRangeByRank(transcript, 1, 5);
        assertNotNull(allCodons);
        GenomeRegion allCodonsPart1 = allCodons.get(0);
        assertEquals(2, allCodonsPart1.start());
        assertEquals(5, allCodonsPart1.end());

        GenomeRegion allCodonsPart2 = allCodons.get(1);
        assertEquals(7, allCodonsPart2.start());
        assertEquals(16, allCodonsPart2.end());

        GenomeRegion allCodonsPart3 = allCodons.get(2);
        assertEquals(18, allCodonsPart3.start());
        assertEquals(18, allCodonsPart3.end());

        assertNull(codonRangeByRank(transcript, 0, 6));
    }

    @Test
    public void canRetrieveCodonRangeByIndexOnReverseGene() {
        HmfTranscriptRegion transcript = create(Strand.REVERSE);

        assertNull(codonRangeByRank(transcript, -1, 1));
        assertNull(codonRangeByRank(transcript, 0, 0));

        List<GenomeRegion> codonsOneAndTwo = codonRangeByRank(transcript, 1, 2);
        assertNotNull(codonsOneAndTwo);
        assertEquals(2, codonsOneAndTwo.size());

        GenomeRegion codonsOneAndTwoPart1 = codonsOneAndTwo.get(0);
        assertEquals(12, codonsOneAndTwoPart1.start());
        assertEquals(16, codonsOneAndTwoPart1.end());

        GenomeRegion codonsOneAndTwoPart2 = codonsOneAndTwo.get(1);
        assertEquals(18, codonsOneAndTwoPart2.start());
        assertEquals(18, codonsOneAndTwoPart2.end());

        List<GenomeRegion> allCodons = codonRangeByRank(transcript, 1, 5);
        assertNotNull(allCodons);
        GenomeRegion allCodonsPart1 = allCodons.get(0);
        assertEquals(2, allCodonsPart1.start());
        assertEquals(5, allCodonsPart1.end());

        GenomeRegion allCodonsPart2 = allCodons.get(1);
        assertEquals(7, allCodonsPart2.start());
        assertEquals(16, allCodonsPart2.end());

        GenomeRegion allCodonsPart3 = allCodons.get(2);
        assertEquals(18, allCodonsPart3.start());
        assertEquals(18, allCodonsPart3.end());

        assertNull(codonRangeByRank(transcript, 0, 6));
    }

    @NotNull
    private static GenomeRegion assertedCodonGet(@NotNull HmfTranscriptRegion transcript, int codonRank) {
        return assertedCodonGet(transcript, codonRank, 0);
    }

    @NotNull
    private static GenomeRegion assertedCodonGet(@NotNull HmfTranscriptRegion transcript, int codonRank, int regionIndex) {
        List<GenomeRegion> regions = codonByIndex(transcript, codonRank);
        assertNotNull(regions);
        return regions.get(regionIndex);
    }

    @NotNull
    private static HmfTranscriptRegion create(@NotNull Strand strand) {
        List<HmfExonRegion> exome = Lists.newArrayList(
                ImmutableHmfExonRegion.builder().chromosome("1").exonRank(1).start(1).end(5).build(),
                ImmutableHmfExonRegion.builder().chromosome("1").exonRank(2).start(7).end(16).build(),
                ImmutableHmfExonRegion.builder().chromosome("1").exonRank(3).start(18).end(20).build());

        return ImmutableHmfTranscriptRegion.builder()
                .chromosome("1")
                .chromosomeBand("band")
                .geneName("gene")
                .geneId("geneID")
                .transName("transcriptID")
                .isCanonical(true)
                .geneStart(1)
                .geneEnd(29)
                .exons(exome)
                .start(1)
                .end(20)
                .codingStart(2)
                .codingEnd(18)
                .strand(strand)
                .build();
    }
}