package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.isofox.ReadCountsTest.REF_BASE_STR_1;
import static com.hartwig.hmftools.isofox.TestUtils.CHR_1;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_1;
import static com.hartwig.hmftools.isofox.TestUtils.createCigar;
import static com.hartwig.hmftools.isofox.TestUtils.createReadPair;
import static com.hartwig.hmftools.isofox.TestUtils.createReadRecord;
import static com.hartwig.hmftools.isofox.TestUtils.createSupplementaryReadPair;
import static com.hartwig.hmftools.isofox.TestUtils.generateRandomBases;
import static com.hartwig.hmftools.isofox.unmapped.UnmappedRead.UMR_NO_MATE;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.GeneTestUtils;
import com.hartwig.hmftools.isofox.common.FragmentTracker;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunctionFinder;
import com.hartwig.hmftools.isofox.unmapped.UmrFinder;
import com.hartwig.hmftools.isofox.unmapped.UnmappedRead;

import org.junit.Test;

import htsjdk.samtools.SAMFlag;

public class UnmappedReadsTest
{
    private UmrFinder mUmrFinder;
    private GeneCollection mGenes;

    public UnmappedReadsTest()
    {
        String chromosome = CHR_1;
        String geneId = GENE_ID_1;

        GeneData geneData = new GeneData(geneId, geneId, chromosome, (byte) 1, 100, 1500, "");

        TranscriptData transData = createTransExons(
                geneId, TRANS_ID_1, POS_STRAND, new int[] {100, 300, 500}, 100, 150, 550, true, "");

        GeneReadData gene = new GeneReadData(geneData);

        gene.setTranscripts(Lists.newArrayList(transData));

        IsofoxConfig config = new IsofoxConfig();
        config.Functions.add(IsofoxFunction.UNMAPPED_READS);
        mUmrFinder = new UmrFinder(config, null);
        mGenes = new GeneCollection(0, Lists.newArrayList(gene));

        mUmrFinder.setGeneData(mGenes);
    }

    @Test
    public void testValidReads()
    {
        //un-paired read
        FragmentTracker fragTracker = new FragmentTracker();

        String refBases = generateRandomBases(70);
        ReadRecord read = createReadRecord(1, CHR_1, 151, 200, refBases, createCigar(0, 50, 20));
        read.setFlag(SAMFlag.MATE_UNMAPPED, true);

        fragTracker.checkRead(read);

        mUmrFinder.processUnpairedReads(fragTracker);

        assertEquals(1, mUmrFinder.getCandidateReads().size());
        UnmappedRead umRead = mUmrFinder.getCandidateReads().get(0);
        assertEquals(SE_END, umRead.ScSide);
        assertEquals(20, umRead.ScLength);
        assertEquals(0, umRead.ExonDistance);
        assertEquals(1, umRead.ExonRank);
        assertTrue(umRead.MateCoords.equals(UMR_NO_MATE));

        mUmrFinder.setGeneData(mGenes); // clear cached data

        // paired reads
        ReadRecord[] reads = createReadPair(1, mGenes, mGenes, 151, 200, 100, 169,
                createCigar(0, 50, 20), createCigar(0, 70, 0), POS_ORIENT, NEG_ORIENT);

        mUmrFinder.processReads(reads[0], reads[1], false);

        assertEquals(1, mUmrFinder.getCandidateReads().size());
        umRead = mUmrFinder.getCandidateReads().get(0);
        assertEquals(SE_END, umRead.ScSide);
        assertEquals(20, umRead.ScLength);
        assertEquals(0, umRead.ExonDistance);
        assertEquals(1, umRead.ExonRank);

        // mark if matching a supplementary alignment
        ReadRecord[] suppReads = createSupplementaryReadPair(1, mGenes, mGenes, 151, 200, 161, 200,
                createCigar(0, 50, 20), createCigar(0, 40, 30), true);

        mUmrFinder.processReads(suppReads[0], suppReads[1], true);
        mUmrFinder.markChimericMatches();

        assertEquals(1, mUmrFinder.getCandidateReads().size());
        assertEquals(1, mUmrFinder.getSupplementaryReadKeys().size());
        assertTrue(umRead.MatchesChimeric);
    }

    @Test
    public void testInvalidReads()
    {
        // reads without soft clips or not sufficient
        ReadRecord[] reads = createReadPair(1, mGenes, mGenes, 120, 191, 300, 359,
                createCigar(0, 70, 0), createCigar(10, 60, 0), POS_ORIENT, NEG_ORIENT);

        mUmrFinder.processReads(reads[0], reads[1], false);

        assertTrue(mUmrFinder.getCandidateReads().isEmpty());

        // soft-clips on both sides of an exon
        reads = createReadPair(1, mGenes, mGenes, 381, 400, 300, 339,
                createCigar(0, 20, 50), createCigar(30, 40, 0), POS_ORIENT, NEG_ORIENT);

        mUmrFinder.processReads(reads[0], reads[1], false);

        assertTrue(mUmrFinder.getCandidateReads().isEmpty());

        // read past the SC boundary
        reads = createReadPair(1, mGenes, mGenes, 381, 400, 450, 519,
                createCigar(0, 20, 50), createCigar(0, 70, 0), POS_ORIENT, NEG_ORIENT);

        mUmrFinder.processReads(reads[0], reads[1], false);

        assertTrue(mUmrFinder.getCandidateReads().isEmpty());

        reads = createReadPair(1, mGenes, mGenes, 250, 319, 300, 339,
                createCigar(0, 70, 50), createCigar(30, 40, 0), POS_ORIENT, NEG_ORIENT);

        mUmrFinder.processReads(reads[0], reads[1], false);

        assertTrue(mUmrFinder.getCandidateReads().isEmpty());

        // more than 5 bases into the intron
        reads = createReadPair(1, mGenes, mGenes, 387, 406, 320, 389,
                createCigar(0, 20, 50), createCigar(0, 70, 0), POS_ORIENT, NEG_ORIENT);

        mUmrFinder.processReads(reads[0], reads[1], false);

        assertTrue(mUmrFinder.getCandidateReads().isEmpty());
    }

}
