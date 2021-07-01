package com.hartwig.hmftools.common.genome.region;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class HmfTranscriptRegionTest {

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

        GenomeRegion codon5Part1 = assertedCodonGet(transcript, 5, 0);
        assertNotNull(codon5Part1);
        assertEquals(15, codon5Part1.start());
        assertEquals(16, codon5Part1.end());

        GenomeRegion codon5Part2 = assertedCodonGet(transcript, 5, 1);
        assertNotNull(codon5Part2);
        assertEquals(18, codon5Part2.start());
        assertEquals(18, codon5Part2.end());

        assertNull(transcript.codonByIndex(6));
    }

    @Test
    public void canRetrieveCodonByIndexOnReverseGene() {
        HmfTranscriptRegion transcript = create(Strand.REVERSE);

        assertNull(transcript.codonByIndex(-1));
        assertNull(transcript.codonByIndex(0));

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
        assertEquals(2, codonsOneAndTwoPart1.start());
        assertEquals(5, codonsOneAndTwoPart1.end());

        GenomeRegion codonsOneAndTwoPart2 = codonsOneAndTwo.get(1);
        assertEquals(7, codonsOneAndTwoPart2.start());
        assertEquals(8, codonsOneAndTwoPart2.end());

        List<GenomeRegion> allCodons = transcript.codonRangeByIndex(1, 5);
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

        assertNull(transcript.codonRangeByIndex(0, 6));
    }

    @Test
    public void canRetrieveCodonRangeByIndexOnReverseGene() {
        HmfTranscriptRegion transcript = create(Strand.REVERSE);

        assertNull(transcript.codonRangeByIndex(-1, 1));
        assertNull(transcript.codonRangeByIndex(0, 0));

        List<GenomeRegion> codonsOneAndTwo = transcript.codonRangeByIndex(1, 2);
        assertNotNull(codonsOneAndTwo);
        assertEquals(2, codonsOneAndTwo.size());

        GenomeRegion codonsOneAndTwoPart1 = codonsOneAndTwo.get(0);
        assertEquals(12, codonsOneAndTwoPart1.start());
        assertEquals(16, codonsOneAndTwoPart1.end());

        GenomeRegion codonsOneAndTwoPart2 = codonsOneAndTwo.get(1);
        assertEquals(18, codonsOneAndTwoPart2.start());
        assertEquals(18, codonsOneAndTwoPart2.end());

        List<GenomeRegion> allCodons = transcript.codonRangeByIndex(1, 5);
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

        assertNull(transcript.codonRangeByIndex(0, 6));
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
        HmfTranscriptRegion braf = HmfGenePanelSupplier.allGenesMap37().get("BRAF");

        List<GenomeRegion> codon600 = braf.codonByIndex(600);
        assertNotNull(codon600);
        assertEquals(1, codon600.size());

        assertEquals(140453135, codon600.get(0).start());
        assertEquals(140453137, codon600.get(0).end());
    }

    @Test
    public void testCodonRangeAtGenomicPosition() {
        HmfTranscriptRegion region = create(Strand.FORWARD);
        assertEquals(0, region.codonRangeAtGenomicPosition(1).size());
        assertEquals(region.codonByIndex(1), region.codonRangeAtGenomicPosition(2));
        assertEquals(region.codonByIndex(1), region.codonRangeAtGenomicPosition(3));
        assertEquals(region.codonByIndex(1), region.codonRangeAtGenomicPosition(4));
        assertEquals(region.codonByIndex(2), region.codonRangeAtGenomicPosition(5));
        assertEquals(0, region.codonRangeAtGenomicPosition(6).size());
        assertEquals(region.codonByIndex(2), region.codonRangeAtGenomicPosition(7));
        assertEquals(region.codonByIndex(2), region.codonRangeAtGenomicPosition(8));
        assertEquals(region.codonByIndex(3), region.codonRangeAtGenomicPosition(9));
        assertEquals(region.codonByIndex(3), region.codonRangeAtGenomicPosition(10));
        assertEquals(region.codonByIndex(3), region.codonRangeAtGenomicPosition(11));
        assertEquals(region.codonByIndex(4), region.codonRangeAtGenomicPosition(12));
        assertEquals(region.codonByIndex(4), region.codonRangeAtGenomicPosition(13));
        assertEquals(region.codonByIndex(4), region.codonRangeAtGenomicPosition(14));
        assertEquals(region.codonByIndex(5), region.codonRangeAtGenomicPosition(15));
        assertEquals(region.codonByIndex(5), region.codonRangeAtGenomicPosition(16));
        assertEquals(0, region.codonRangeAtGenomicPosition(17).size());
        assertEquals(region.codonByIndex(5), region.codonRangeAtGenomicPosition(18));
        assertEquals(0, region.codonRangeAtGenomicPosition(19).size());
        assertEquals(0, region.codonRangeAtGenomicPosition(20).size());
        assertEquals(0, region.codonRangeAtGenomicPosition(21).size());
    }

    @Test
    public void testCodonRangeAtGenomicPositionWorksForRealNRAS() {
        HmfTranscriptRegion nras = HmfGenePanelSupplier.allGenesMap37().get("NRAS");
        assertEquals(0, nras.codonRangeAtGenomicPosition(115251155).size());
        assertEquals(nras.codonByIndex(190), nras.codonRangeAtGenomicPosition(115251156));
        assertEquals(nras.codonByIndex(190), nras.codonRangeAtGenomicPosition(115251157));
        assertEquals(nras.codonByIndex(190), nras.codonRangeAtGenomicPosition(115251158));
        assertEquals(nras.codonByIndex(189), nras.codonRangeAtGenomicPosition(115251159));
        assertEquals(nras.codonByIndex(189), nras.codonRangeAtGenomicPosition(115251160));
        assertEquals(nras.codonByIndex(189), nras.codonRangeAtGenomicPosition(115251161));
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
        List<HmfExonRegion> exome = Lists.newArrayList(
                ImmutableHmfExonRegion.builder().chromosome("1").exonRank(1).start(1).end(5).build(),
                ImmutableHmfExonRegion.builder().chromosome("1").exonRank(2).start(7).end(16).build(),
                ImmutableHmfExonRegion.builder().chromosome("1").exonRank(3).start(18).end(20).build());

        return ImmutableHmfTranscriptRegion.builder()
                .chromosome("1")
                .chromosomeBand("band")
                .gene("gene")
                .geneID("geneID")
                .transcriptID("transcriptID")
                .transcriptVersion(1)
                .geneStart(1)
                .geneEnd(29)
                .exome(exome)
                .start(1)
                .end(20)
                .codingStart(2)
                .codingEnd(18)
                .strand(strand)
                .build();
    }
}