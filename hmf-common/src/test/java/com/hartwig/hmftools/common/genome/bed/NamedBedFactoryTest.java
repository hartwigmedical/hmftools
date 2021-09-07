package com.hartwig.hmftools.common.genome.bed;

import static junit.framework.TestCase.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.ImmutableHmfExonRegion;
import com.hartwig.hmftools.common.genome.region.ImmutableHmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.Strand;

import org.junit.Test;

public class NamedBedFactoryTest {

    private static final String CONTIG = "1";

    @Test
    public void testForwardCodingRegionsWithoutUTR() {
        List<NamedBed> bed = NamedBedFactory.codingRegions(false, create());
        assertEquals(3, bed.size());
        assertNamedBed(bed.get(0), 1050, 1110);
        assertNamedBed(bed.get(1), 1190, 1310);
        assertNamedBed(bed.get(2), 1390, 1450);
    }

    @Test
    public void testForwardCodingRegionsWithUTR() {
        List<NamedBed> bed = NamedBedFactory.codingRegions(true, create());
        assertEquals(3, bed.size());
        assertNamedBed(bed.get(0), 1000, 1110);
        assertNamedBed(bed.get(1), 1190, 1310);
        assertNamedBed(bed.get(2), 1390, 1500);
    }

    private void assertNamedBed(NamedBed victim, long expectedStart, long expectedEnd) {
        assertEquals(expectedStart, victim.start());
        assertEquals(expectedEnd, victim.end());
    }

    private HmfTranscriptRegion create() {

        HmfExonRegion first = exon(1000, 1100);
        HmfExonRegion second = exon(1200, 1300);
        HmfExonRegion third = exon(1400, 1500);

        return ImmutableHmfTranscriptRegion.builder()
                .chromosome(CONTIG)
                .start(1000)
                .end(1500)
                .chromosomeBand("band")
                .gene("gene")
                .geneStart(1000)
                .geneEnd(1500)
                .codingStart(1050)
                .codingEnd(1450)
                .addExome(first)
                .addExome(second)
                .addExome(third)
                .strand(Strand.REVERSE)
                .geneID("geneId")
                .transcriptID("transcriptId")
                .build();

    }

    private HmfExonRegion exon(int start, int end) {
        return ImmutableHmfExonRegion.builder().chromosome(CONTIG).start(start).end(end).exonRank(1).build();
    }

}
