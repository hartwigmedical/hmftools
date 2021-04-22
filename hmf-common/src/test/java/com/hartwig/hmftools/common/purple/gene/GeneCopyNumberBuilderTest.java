package com.hartwig.hmftools.common.purple.gene;

import static org.junit.Assert.assertEquals;

import java.util.Collections;

import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.ImmutableHmfExonRegion;
import com.hartwig.hmftools.common.genome.region.ImmutableHmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.purple.PurpleTestUtils;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

public class GeneCopyNumberBuilderTest {

    private static final String CHROMOSOME = "1";
    private static final double EPSILON = 1E-10;

    private GeneCopyNumberBuilder victim;

    @Before
    public void setup() {
        HmfTranscriptRegion gene = create(1001, 10000);
        victim = new GeneCopyNumberBuilder(gene);
    }

    @Test
    public void testOneCopyNumberOneExon() {
        addCopyNumber(1, 10000, 2);
        addExon(1001, 2000);
        assertCopyNumber(1, 1, 2, 2);
    }

    @Test
    public void testOneCopyNumberTwoExons() {
        addCopyNumber(1, 10000, 2);
        addExon(1001, 2000);
        addExon(3001, 4000);
        assertCopyNumber(1, 1, 2, 2);
    }

    @Test
    public void testAverageOfExonBases() {
        addCopyNumber(1, 1500, 2);
        addExon(1001, 2000);
        addCopyNumber(1501, 10000, 3);
        addExon(3001, 4000);
        assertCopyNumber(2, 1, 2, 3);
    }

    @Test
    public void testCopyNumberChangeInIntron() {
        addCopyNumber(1, 2500, 2);
        addExon(1001, 2000);
        addCopyNumber(2501, 3000, 3);
        addCopyNumber(3001, 10000, 2);
        addExon(3001, 4000);
        assertCopyNumber(1, 1, 2, 2);
    }

    @Test
    public void testGermlineAmplificationInExon() {
        addCopyNumber(1, 2500, 2);
        addExon(1001, 3000);
        addCopyNumber(2501, 3000, 3);
        addCopyNumber(3001, 10000, 2);
        addExon(3001, 4000);
        assertCopyNumber(3, 2, 2, 3);
    }

    @Test
    public void testSingleNegativeRegion() {
        addCopyNumber(1, 2500, -0.8);
        addExon(1001, 2000);
        assertCopyNumber(1, 1, -0.8, -0.8);
    }

    @Test
    public void testSingleZeroRegion() {
        addCopyNumber(1, 2500, 0);
        addExon(1001, 2000);
        assertCopyNumber(1, 1, 0, 0);
    }

    private void assertCopyNumber(int somaticCount, final int minCount, double expectedMin, double expectedMax) {
        final GeneCopyNumber geneCopyNumber = victim.build();
        assertEquals(somaticCount, geneCopyNumber.somaticRegions());
        assertEquals(expectedMin, geneCopyNumber.minCopyNumber(), EPSILON);
        assertEquals(expectedMax, geneCopyNumber.maxCopyNumber(), EPSILON);
    }

    private void addExon(long start, long end) {
        victim.secondary(exon(start, end));
    }

    private void addCopyNumber(long start, long end, double copyNumber) {
        victim.primary(createCopyNumber(start, end, copyNumber));
    }

    @NotNull
    private static PurpleCopyNumber createCopyNumber(long start, long end, double copyNumber) {
        return PurpleTestUtils.createCopyNumber(CHROMOSOME, start, end, copyNumber).build();
    }

    @NotNull
    private static HmfExonRegion exon(long start, long end) {
        return ImmutableHmfExonRegion.builder().exonID("ID").chromosome(CHROMOSOME).start(start).end(end).build();
    }

    @NotNull
    private static HmfTranscriptRegion create(long start, long end) {
        return ImmutableHmfTranscriptRegion.builder()
                .chromosome(CHROMOSOME)
                .start(start)
                .end(end)
                .gene("GENE")
                .transcriptID("ID")
                .transcriptVersion(1)
                .chromosomeBand("BAND")
                .entrezId(Collections.singletonList(1))
                .geneID("ID")
                .geneStart(start)
                .geneEnd(end)
                .codingStart(0)
                .codingEnd(0)
                .strand(Strand.FORWARD)
                .build();
    }
}
