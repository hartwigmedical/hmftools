package com.hartwig.hmftools.common.gene;

import static org.junit.Assert.assertEquals;

import java.util.Collections;

import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.region.hmfslicer.HmfExonRegion;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.region.hmfslicer.ImmutableHmfExonRegion;
import com.hartwig.hmftools.common.region.hmfslicer.ImmutableHmfGenomeRegion;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

public class GeneCopyNumberBuilderTest {

    private static final String CHROMOSOME = "1";
    private static final double EPSILON = 1E-10;

    private HmfGenomeRegion gene;
    private GeneCopyNumberBuilder victim;

    @Before
    public void setup() {
        gene = create(1001, 10000);
        victim = new GeneCopyNumberBuilder(gene);
    }

    @Test
    public void testOneCopyNumberOneExon() {
        addCopyNumber(1, 10000, 2);
        addExon(1001, 2000);
        assertCopyNumber(1, 2, 2, 2);
    }

    @Test
    public void testOneCopyNumberTwoExons() {
        addCopyNumber(1, 10000, 2);
        addExon(1001, 2000);
        addExon(3001, 4000);
        assertCopyNumber(1, 2, 2, 2);
    }

    @Test
    public void testAverageOfExonBases() {
        addCopyNumber(1, 1500, 2);
        addExon(1001, 2000);
        addCopyNumber(1501, 10000, 3);
        addExon(3001, 4000);
        assertCopyNumber(2, 2, 2.75, 3);
    }

    @Test
    public void testCopyNumberChangeInIntron() {
        addCopyNumber(1, 2500, 2);
        addExon(1001, 2000);
        addCopyNumber(2501, 3000, 3);
        addCopyNumber(3001, 10000, 2);
        addExon(3001, 4000);
        assertCopyNumber(1, 2, 2, 2);
    }

    @Test
    public void testSingleNegativeRegion() {
        addCopyNumber(1, 2500, -0.8);
        addExon(1001, 2000);
        assertCopyNumber(1, -0.8, -0.8, -0.8);
    }

    private void assertCopyNumber(int count, double expectedMin, double expectedMean, double expectedMax) {
        final GeneCopyNumber geneCopyNumber = victim.build();
        assertEquals(count, geneCopyNumber.regions());
        assertEquals(expectedMin, geneCopyNumber.minCopyNumber(), EPSILON);
        assertEquals(expectedMean, geneCopyNumber.meanCopyNumber(), EPSILON);
        assertEquals(expectedMax, geneCopyNumber.maxCopyNumber(), EPSILON);
    }

    private void addExon(long start, long end) {
        victim.addExon(exon(start, end));
    }

    private void addCopyNumber(long start, long end, double copyNumber) {
        victim.addCopyNumber(createCopyNumber(start, end, copyNumber));
    }

    @NotNull
    private static PurpleCopyNumber createCopyNumber(long start, long end, double copyNumber) {
        return PurpleDatamodelTest.createCopyNumber(CHROMOSOME, start, end, copyNumber).build();
    }

    @NotNull
    private static HmfExonRegion exon(long start, long end) {
        return ImmutableHmfExonRegion.builder().exonID("ID").chromosome(CHROMOSOME).start(start).end(end).build();
    }

    private static HmfGenomeRegion create(long start, long end) {
        return ImmutableHmfGenomeRegion.builder()
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
                .build();
    }

}
