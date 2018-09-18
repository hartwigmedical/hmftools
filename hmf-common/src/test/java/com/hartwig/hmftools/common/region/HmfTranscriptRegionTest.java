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
    //codons:            |___________|_______________|_____________|______________|

    @Test
    public void canRetrieveExonByIndex() {
        HmfTranscriptRegion forwardStrandGene = create(Strand.FORWARD);
        HmfExonRegion exon0 = forwardStrandGene.exonByIndex(0);
        assertNull(exon0);

        HmfExonRegion exon1 = forwardStrandGene.exonByIndex(1);
        assertNotNull(exon1);
        assertEquals("1", exon1.exonID());

        HmfExonRegion exon2 = forwardStrandGene.exonByIndex(2);
        assertNotNull(exon2);
        assertEquals("2", exon2.exonID());

        HmfExonRegion exon3 = forwardStrandGene.exonByIndex(3);
        assertNull(exon3);

        HmfTranscriptRegion reverseStrandGene = create(Strand.REVERSE);
        exon0 = reverseStrandGene.exonByIndex(0);
        assertNull(exon0);

        exon1 = reverseStrandGene.exonByIndex(1);
        assertNotNull(exon1);
        assertEquals("2", exon1.exonID());

        exon2 = reverseStrandGene.exonByIndex(2);
        assertNotNull(exon2);
        assertEquals("1", exon2.exonID());

        exon3 = reverseStrandGene.exonByIndex(3);
        assertNull(exon3);
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