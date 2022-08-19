package com.hartwig.hmftools.purple.gene;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_2;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.purple.PurpleTestUtils;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.test.GeneTestUtils;

import org.junit.Test;

public class GeneCopyNumberBuilderTest
{
    private static final String CHROMOSOME = "1";
    private static final double EPSILON = 1E-10;

    private final GeneData mGeneData = GeneTestUtils.createEnsemblGeneData(
            GENE_ID_1, GENE_NAME_1, CHR_1, POS_STRAND,1001, 10000);

    private TranscriptData mTransData = GeneTestUtils.createTransExons(
            GENE_ID_1, TRANS_ID_1, POS_STRAND, new int[] { 1001, 3001}, 999, null, null, true, "");

    private TranscriptData mTransDataSingleExon = GeneTestUtils.createTransExons(
            GENE_ID_1, TRANS_ID_2, POS_STRAND, new int[] { 1001 }, 999, null, null, true, "");

    @Test
    public void testOneCopyNumberOneExon()
    {
        List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();
        copyNumbers.add(createCopyNumber(1, 10000, 2));

        GeneCopyNumberBuilder builder = new GeneCopyNumberBuilder(mGeneData, mTransDataSingleExon, copyNumbers);
        GeneCopyNumber geneCopyNumber = builder.create();
        assertCopyNumber(geneCopyNumber, 1,  2, 2);
    }

    @Test
    public void testOneCopyNumberTwoExons()
    {
        List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();
        copyNumbers.add(createCopyNumber(1, 10000, 2));

        GeneCopyNumberBuilder builder = new GeneCopyNumberBuilder(mGeneData, mTransData, copyNumbers);
        GeneCopyNumber geneCopyNumber = builder.create();
        assertCopyNumber(geneCopyNumber, 1, 2, 2);
    }

    @Test
    public void testAverageOfExonBases()
    {
        List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();
        copyNumbers.add(createCopyNumber(1, 1500, 2));
        copyNumbers.add(createCopyNumber(1501, 10000, 3));

        GeneCopyNumberBuilder builder = new GeneCopyNumberBuilder(mGeneData, mTransData, copyNumbers);
        GeneCopyNumber geneCopyNumber = builder.create();
        assertCopyNumber(geneCopyNumber, 2, 2, 3);
    }

    @Test
    public void testCopyNumberChangeInIntron()
    {
        List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();
        copyNumbers.add(createCopyNumber(1, 2500, 2));
        copyNumbers.add(createCopyNumber(2501, 3000, 3));
        copyNumbers.add(createCopyNumber(3001, 10000, 2));

        GeneCopyNumberBuilder builder = new GeneCopyNumberBuilder(mGeneData, mTransDataSingleExon, copyNumbers);
        GeneCopyNumber geneCopyNumber = builder.create();
        assertCopyNumber(geneCopyNumber, 1,  2, 2);
    }

    @Test
    public void testSingleNegativeRegion()
    {
        List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();
        copyNumbers.add(createCopyNumber(1, 2500, -0.8));

        GeneCopyNumberBuilder builder = new GeneCopyNumberBuilder(mGeneData, mTransDataSingleExon, copyNumbers);
        GeneCopyNumber geneCopyNumber = builder.create();
        assertCopyNumber(geneCopyNumber, 1, -0.8, -0.8);
    }

    @Test
    public void testSingleZeroRegion()
    {
        List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();
        copyNumbers.add(createCopyNumber(1, 2500, 0));

        GeneCopyNumberBuilder builder = new GeneCopyNumberBuilder(mGeneData, mTransDataSingleExon, copyNumbers);
        GeneCopyNumber geneCopyNumber = builder.create();
        assertCopyNumber(geneCopyNumber, 1, 0, 0);
    }

    private void assertCopyNumber(final GeneCopyNumber geneCopyNumber, int somaticCount, double expectedMin, double expectedMax)
    {
        assertEquals(somaticCount, geneCopyNumber.somaticRegions());
        assertEquals(expectedMin, geneCopyNumber.minCopyNumber(), EPSILON);
        assertEquals(expectedMax, geneCopyNumber.maxCopyNumber(), EPSILON);
    }

    private static PurpleCopyNumber createCopyNumber(int start, int end, double copyNumber)
    {
        return PurpleTestUtils.createCopyNumber(CHROMOSOME, start, end, copyNumber).build();
    }
}
