package com.hartwig.hmftools.common.genome.region;

import static com.hartwig.hmftools.common.genome.region.HmfTranscriptRegionUtils.codonByIndex;
import static com.hartwig.hmftools.common.genome.region.HmfTranscriptRegionUtils.codonRangeAtGenomicPosition;
import static com.hartwig.hmftools.common.genome.region.HmfTranscriptRegionUtils.codonRangeByIndex;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;

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

        assertNull(codonRangeByIndex(transcript, -1, 1));
        assertNull(codonRangeByIndex(transcript, 0, 0));

        List<GenomeRegion> codonsOneAndTwo = codonRangeByIndex(transcript, 1, 2);
        assertNotNull(codonsOneAndTwo);
        GenomeRegion codonsOneAndTwoPart1 = codonsOneAndTwo.get(0);
        assertEquals(2, codonsOneAndTwoPart1.start());
        assertEquals(5, codonsOneAndTwoPart1.end());

        GenomeRegion codonsOneAndTwoPart2 = codonsOneAndTwo.get(1);
        assertEquals(7, codonsOneAndTwoPart2.start());
        assertEquals(8, codonsOneAndTwoPart2.end());

        List<GenomeRegion> allCodons = codonRangeByIndex(transcript, 1, 5);
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

        assertNull(codonRangeByIndex(transcript, 0, 6));
    }

    @Test
    public void canRetrieveCodonRangeByIndexOnReverseGene() {
        HmfTranscriptRegion transcript = create(Strand.REVERSE);

        assertNull(codonRangeByIndex(transcript, -1, 1));
        assertNull(codonRangeByIndex(transcript, 0, 0));

        List<GenomeRegion> codonsOneAndTwo = codonRangeByIndex(transcript, 1, 2);
        assertNotNull(codonsOneAndTwo);
        assertEquals(2, codonsOneAndTwo.size());

        GenomeRegion codonsOneAndTwoPart1 = codonsOneAndTwo.get(0);
        assertEquals(12, codonsOneAndTwoPart1.start());
        assertEquals(16, codonsOneAndTwoPart1.end());

        GenomeRegion codonsOneAndTwoPart2 = codonsOneAndTwo.get(1);
        assertEquals(18, codonsOneAndTwoPart2.start());
        assertEquals(18, codonsOneAndTwoPart2.end());

        List<GenomeRegion> allCodons = codonRangeByIndex(transcript, 1, 5);
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

        assertNull(codonRangeByIndex(transcript, 0, 6));
    }

    @Test
    public void worksForRealBRAFCodon600() {
        HmfTranscriptRegion braf = HmfGenePanelSupplier.allGenesMap37().get("BRAF");

        List<GenomeRegion> codon600 = codonByIndex(braf, 600);
        assertNotNull(codon600);
        assertEquals(1, codon600.size());

        assertEquals(140453135, codon600.get(0).start());
        assertEquals(140453137, codon600.get(0).end());
    }

    @Test
    public void testCodonRangeAtGenomicPosition() {
        HmfTranscriptRegion region = create(Strand.FORWARD);
        assertEquals(0, codonRangeAtGenomicPosition(region, 1).size());
        assertEquals(codonByIndex(region, 1), codonRangeAtGenomicPosition(region, 2));
        assertEquals(codonByIndex(region, 1), codonRangeAtGenomicPosition(region, 3));
        assertEquals(codonByIndex(region, 1), codonRangeAtGenomicPosition(region, 4));
        assertEquals(codonByIndex(region, 2), codonRangeAtGenomicPosition(region, 5));
        assertEquals(0, codonRangeAtGenomicPosition(region, 6).size());
        assertEquals(codonByIndex(region, 2), codonRangeAtGenomicPosition(region, 7));
        assertEquals(codonByIndex(region, 2), codonRangeAtGenomicPosition(region, 8));
        assertEquals(codonByIndex(region, 3), codonRangeAtGenomicPosition(region, 9));
        assertEquals(codonByIndex(region, 3), codonRangeAtGenomicPosition(region, 10));
        assertEquals(codonByIndex(region, 3), codonRangeAtGenomicPosition(region, 11));
        assertEquals(codonByIndex(region, 4), codonRangeAtGenomicPosition(region, 12));
        assertEquals(codonByIndex(region, 4), codonRangeAtGenomicPosition(region, 13));
        assertEquals(codonByIndex(region, 4), codonRangeAtGenomicPosition(region, 14));
        assertEquals(codonByIndex(region, 5), codonRangeAtGenomicPosition(region, 15));
        assertEquals(codonByIndex(region, 5), codonRangeAtGenomicPosition(region, 16));
        assertEquals(0, codonRangeAtGenomicPosition(region, 17).size());
        assertEquals(codonByIndex(region, 5), codonRangeAtGenomicPosition(region, 18));
        assertEquals(0, codonRangeAtGenomicPosition(region, 19).size());
        assertEquals(0, codonRangeAtGenomicPosition(region, 20).size());
        assertEquals(0, codonRangeAtGenomicPosition(region, 21).size());
    }

    @Test
    public void testCodonRangeAtGenomicPositionWorksForRealNRAS() {
        HmfTranscriptRegion nras = HmfGenePanelSupplier.allGenesMap37().get("NRAS");
        assertEquals(0, codonRangeAtGenomicPosition(nras, 115251155).size());
        assertEquals(codonByIndex(nras, 190), codonRangeAtGenomicPosition(nras, 115251156));
        assertEquals(codonByIndex(nras, 190), codonRangeAtGenomicPosition(nras, 115251157));
        assertEquals(codonByIndex(nras, 190), codonRangeAtGenomicPosition(nras, 115251158));
        assertEquals(codonByIndex(nras, 189), codonRangeAtGenomicPosition(nras, 115251159));
        assertEquals(codonByIndex(nras, 189), codonRangeAtGenomicPosition(nras, 115251160));
        assertEquals(codonByIndex(nras, 189), codonRangeAtGenomicPosition(nras, 115251161));
    }

    @NotNull
    private static GenomeRegion assertedCodonGet(@NotNull HmfTranscriptRegion transcript, int codonIndex) {
        return assertedCodonGet(transcript, codonIndex, 0);
    }

    @NotNull
    private static GenomeRegion assertedCodonGet(@NotNull HmfTranscriptRegion transcript, int codonIndex, int regionIndex) {
        List<GenomeRegion> regions = codonByIndex(transcript, codonIndex);
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