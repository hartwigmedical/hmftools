package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createGeneDataCache;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.isofox.IsofoxConstants.MAX_NOVEL_SJ_DISTANCE;
import static com.hartwig.hmftools.isofox.TestUtils.CHR_1;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_1;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_2;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_5;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_6;
import static com.hartwig.hmftools.isofox.TestUtils.addTestGenes;
import static com.hartwig.hmftools.isofox.TestUtils.addTestTranscripts;
import static com.hartwig.hmftools.isofox.TestUtils.createCigar;
import static com.hartwig.hmftools.isofox.TestUtils.createGeneCollection;
import static com.hartwig.hmftools.isofox.TestUtils.createMappedRead;
import static com.hartwig.hmftools.isofox.TestUtils.createSupplementaryReadPair;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.MATCHED_JUNCTION;

import static org.junit.Assert.assertEquals;

import static htsjdk.samtools.SAMFlag.FIRST_OF_PAIR;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.isofox.common.BaseDepth;
import com.hartwig.hmftools.isofox.common.FragmentTracker;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.fusion.ChimericReadTracker;

import org.junit.Test;

public class ChimericReadTest
{
    /*
    Single-gene tests
        - INV, BND straight to chimeric
        - DEL - chimeric if 500K up/down stream
        - DEL if between diff genes by known splice donor/acceptor


    */

    @Test
    public void testBasicReads()
    {
        final EnsemblDataCache geneTransCache = createGeneDataCache();

        addTestGenes(geneTransCache);
        addTestTranscripts(geneTransCache);

        int gcId = 0;

        final GeneCollection gc1 =
                createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_1)));

        IsofoxConfig config = new IsofoxConfig();
        ChimericReadTracker chimericRT = new ChimericReadTracker(config);

        chimericRT.initialise(gc1);
        BaseDepth baseDepth = new BaseDepth();
        baseDepth.initialise(gc1.regionBounds());

        FragmentTracker fragTracker = new FragmentTracker();

        // chimeric read pair (eg a BND)
        int readId = 0;
        ReadRecord read1 = createMappedRead(readId, gc1, 1050, 1089, createCigar(0, 40, 0));
        read1.setFlag(FIRST_OF_PAIR, true);
        ReadRecord read2 = createMappedRead(readId, gc1, 1081, 1100, createCigar(0, 20, 20));
        read2.setSuppAlignment("supp");

        chimericRT.addChimericReadPair(read1, read2);
        chimericRT.postProcessChimericReads(baseDepth, fragTracker);

        assertEquals(1, chimericRT.getReadMap().size());
        assertEquals(2, chimericRT.getReadMap().get(read1.Id).size());
        assertEquals(1, chimericRT.getJunctionPositions().size());
        assertEquals(1100, chimericRT.getJunctionPositions().iterator().next().intValue());

        // single read in this gene
        chimericRT.clear();

        ReadRecord read3 = createMappedRead(++readId, gc1, 1081, 1100, createCigar(0, 20, 20));
        read3.setSuppAlignment("supp");
        fragTracker.checkRead(read3);

        // realignable read pair supporting a junction
        read1 = createMappedRead(++readId, gc1, 1050, 1089, createCigar(0, 40, 0));
        read1.setFlag(FIRST_OF_PAIR, true);
        read2 = createMappedRead(readId, gc1, 1066, 1100, createCigar(0, 35, 5));

        chimericRT.addRealignmentCandidates(read1, read2);

        chimericRT.postProcessChimericReads(baseDepth, fragTracker);

        assertEquals(2, chimericRT.getReadMap().size());
        assertEquals(1, chimericRT.getReadMap().get(read3.Id).size());
        assertEquals(2, chimericRT.getReadMap().get(read1.Id).size());
        assertEquals(1, chimericRT.getJunctionPositions().size());
        assertEquals(1100, chimericRT.getJunctionPositions().iterator().next().intValue());
    }

    @Test
    public void testSameGeneCollection()
    {
        final EnsemblDataCache geneTransCache = createGeneDataCache();

        addTestGenes(geneTransCache);
        addTestTranscripts(geneTransCache);

        int gcId = 0;

        // test with 2 overlapping genes
        List<EnsemblGeneData> geneDataList = Lists.newArrayList(
                geneTransCache.getGeneDataById(GENE_ID_1), geneTransCache.getGeneDataById(GENE_ID_2));

        final GeneCollection gc1 = createGeneCollection(geneTransCache, gcId++, geneDataList);

        IsofoxConfig config = new IsofoxConfig();
        ChimericReadTracker chimericRT = new ChimericReadTracker(config);

        chimericRT.initialise(gc1);
        BaseDepth baseDepth = new BaseDepth();
        baseDepth.initialise(gc1.regionBounds());

        FragmentTracker fragTracker = new FragmentTracker();

        // DEL only covering 1 gene isn't considered chimeric
        int readId = 0;
        ReadRecord read1 = createMappedRead(readId, gc1, 1081, 1100, createCigar(0, 20, 20));
        read1.setFlag(FIRST_OF_PAIR, true);
        ReadRecord read2 = createMappedRead(readId, gc1, 1500, 1519, createCigar(20, 20, 0));
        read2.setStrand(true, false);

        chimericRT.addChimericReadPair(read1, read2);
        chimericRT.postProcessChimericReads(baseDepth, fragTracker);

        assertTrue(chimericRT.getReadMap().isEmpty());
        assertEquals(1, chimericRT.getLocalChimericReads().size());
        assertTrue(chimericRT.getJunctionPositions().isEmpty());

        chimericRT.clear();

        // DEL linking 2 genes at known splice sites is considered chimeric
        read1 = createMappedRead(++readId, gc1, 1081, 1100, createCigar(0, 20, 20));
        read1.setFlag(FIRST_OF_PAIR, true);
        read2 = createMappedRead(readId, gc1, 10500, 10519, createCigar(20, 20, 0));
        read2.setStrand(true, false);

        chimericRT.addChimericReadPair(read1, read2);
        chimericRT.postProcessChimericReads(baseDepth, fragTracker);

        assertEquals(1, chimericRT.getReadMap().size());
        assertEquals(2, chimericRT.getJunctionPositions().size());
        assertTrue(chimericRT.getJunctionPositions().contains(1100));
        assertTrue(chimericRT.getJunctionPositions().contains(10500));

        chimericRT.clear();

        // this DEL doesn't splice at known sites
        read1 = createMappedRead(++readId, gc1, 1081, 1100, createCigar(0, 20, 20));
        read1.setFlag(FIRST_OF_PAIR, true);
        read2 = createMappedRead(readId, gc1, 10450, 10469, createCigar(20, 20, 0));
        read2.setStrand(true, false);

        chimericRT.addChimericReadPair(read1, read2);
        chimericRT.postProcessChimericReads(baseDepth, fragTracker);

        assertTrue(chimericRT.getReadMap().isEmpty());
        assertEquals(1, chimericRT.getLocalChimericReads().size());
        assertTrue(chimericRT.getJunctionPositions().isEmpty());

        chimericRT.clear();

        // unless the junction could be explained by a single gene, as is the case for 10500 below
        geneDataList = Lists.newArrayList(
                geneTransCache.getGeneDataById(GENE_ID_5), geneTransCache.getGeneDataById(GENE_ID_6));

        final GeneCollection gc2 = createGeneCollection(geneTransCache, gcId++, geneDataList);

        read1 = createMappedRead(++readId, gc2, 10481, 10500, createCigar(0, 20, 20));
        read1.setFlag(FIRST_OF_PAIR, true);
        read2 = createMappedRead(readId, gc2, 11200, 11219, createCigar(20, 20, 0));
        read2.setStrand(true, false);

        chimericRT.initialise(gc2);

        chimericRT.addChimericReadPair(read1, read2);
        chimericRT.postProcessChimericReads(baseDepth, fragTracker);

        assertTrue(chimericRT.getReadMap().isEmpty());
        assertEquals(1, chimericRT.getLocalChimericReads().size());
        assertTrue(chimericRT.getJunctionPositions().isEmpty());

        chimericRT.clear();

        // a DEL more than 500K from the gene boundary is considered chimeric
        read1 = createMappedRead(++readId, gc2, 10481, 10500, createCigar(0, 20, 20));
        read1.setFlag(FIRST_OF_PAIR, true);
        int geneEnd = gc2.genes().stream().mapToInt(x -> x.GeneData.GeneEnd).max().orElse(0);
        int chimericJunc = geneEnd + MAX_NOVEL_SJ_DISTANCE + 1;
        read2 = createMappedRead(readId, gc2, chimericJunc, chimericJunc + 19, createCigar(20, 20, 0));
        read2.setStrand(true, false);
        read2.setGeneCollection(SE_START, gc2.id(), false);
        read2.setGeneCollection(SE_END, gc2.id(), false);

        chimericRT.addChimericReadPair(read1, read2);
        chimericRT.postProcessChimericReads(baseDepth, fragTracker);

        assertEquals(1, chimericRT.getReadMap().size());
        assertEquals(2, chimericRT.getJunctionPositions().size());
        assertTrue(chimericRT.getJunctionPositions().contains(10500));
        assertTrue(chimericRT.getJunctionPositions().contains(chimericJunc));
    }

    @Test
    public void testPrePosGeneReads()
    {
        // 2 gene collections, testing reads before, inside and after the genes
        final EnsemblDataCache geneTransCache = createGeneDataCache();

        addTestGenes(geneTransCache);
        addTestTranscripts(geneTransCache);

        int gcId = 0;

        final GeneCollection gc1 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_1)));
        final GeneCollection gc2 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_2)));

        IsofoxConfig config = new IsofoxConfig();
        ChimericReadTracker chimericRT = new ChimericReadTracker(config);
        BaseDepth baseDepth = new BaseDepth();

        chimericRT.initialise(gc1);
        baseDepth.initialise(gc1.regionBounds());

        FragmentTracker fragTracker = new FragmentTracker();

        // pre and post reads are discarded since they don't affect the gene
        int readId = 0;
        ReadRecord read1 = createMappedRead(readId, gc1, 481, 500, createCigar(0, 20, 20));
        read1.setFlag(FIRST_OF_PAIR, true);
        ReadRecord read2 = createMappedRead(readId, gc1, 2000, 2019, createCigar(20, 20, 0));
        read2.setStrand(true, false);

        chimericRT.addChimericReadPair(read1, read2);
        chimericRT.postProcessChimericReads(baseDepth, fragTracker);

        assertTrue(chimericRT.getReadMap().isEmpty());
        assertTrue(chimericRT.getLocalChimericReads().isEmpty());
        assertTrue(chimericRT.getJunctionPositions().isEmpty());
        chimericRT.clear();

        chimericRT.initialise(gc2);
        baseDepth.initialise(gc2.regionBounds());
        fragTracker.checkRead(read2);
        chimericRT.postProcessChimericReads(baseDepth, fragTracker);

        assertTrue(chimericRT.getReadMap().isEmpty());
        assertTrue(chimericRT.getLocalChimericReads().isEmpty());
        assertTrue(chimericRT.getJunctionPositions().isEmpty());

        // pre and post gene reads are kept if relating to other genes / have supp alignments
        chimericRT.clear();
        fragTracker.clear();
        chimericRT.initialise(gc1);
        baseDepth.initialise(gc1.regionBounds());

        read1 = createMappedRead(++readId, gc1, 481, 500, createCigar(0, 20, 20));
        read1.setFlag(FIRST_OF_PAIR, true);
        read1.setSuppAlignment("supp");

        read2 = createMappedRead(++readId, gc1, 2000, 2019, createCigar(20, 20, 0));
        read2.setStrand(true, false);
        read2.setSuppAlignment("supp");

        // these post-gene reads will be skipped in this gene collection
        ReadRecord read3 = createMappedRead(readId, gc1, 2010, 2049, createCigar(0, 40, 0));

        // will also be skipped since relates to next gene collection
        ReadRecord read4 = createMappedRead(++readId, gc1, 2000, 2019, createCigar(20, 20, 0));

        // will be processed (as local alt-SJ candidates) and the post gene read cached so as not to handle again in gc2
        ReadRecord read5 = createMappedRead(++readId, gc1, 1081, 1100, createCigar(0, 20, 20));
        ReadRecord read6 = createMappedRead(readId, gc1, 2100, 2119, createCigar(20, 20, 0));

        fragTracker.checkRead(read1);
        fragTracker.checkRead(read4);
        chimericRT.addChimericReadPair(read2, read3);
        chimericRT.addChimericReadPair(read5, read6);

        chimericRT.postProcessChimericReads(baseDepth, fragTracker);

        assertEquals(1, chimericRT.getReadMap().size());
        assertTrue(chimericRT.getReadMap().containsKey(read1.Id));
        assertFalse(chimericRT.getReadMap().containsKey(read2.Id));
        assertFalse(chimericRT.getReadMap().containsKey(read4.Id));
        assertFalse(chimericRT.getReadMap().containsKey(read5.Id));

        assertEquals(1, chimericRT.getLocalChimericReads().size());
        assertTrue(chimericRT.getLocalChimericReads().get(0).contains(read5));

        assertTrue(chimericRT.getJunctionPositions().contains(500));
        assertFalse(chimericRT.getJunctionPositions().contains(2000));

        chimericRT.clear();
        fragTracker.clear();
        chimericRT.initialise(gc2);
        baseDepth.initialise(gc2.regionBounds());

        chimericRT.addChimericReadPair(read2, read3);
        fragTracker.checkRead(read4);
        fragTracker.checkRead(read6);

        chimericRT.postProcessChimericReads(baseDepth, fragTracker);

        assertEquals(2, chimericRT.getReadMap().size());
        assertTrue(chimericRT.getReadMap().containsKey(read2.Id));
        assertTrue(chimericRT.getReadMap().containsKey(read4.Id));
        assertTrue(chimericRT.getLocalChimericReads().isEmpty());
        assertTrue(chimericRT.getJunctionPositions().contains(2000));
        assertFalse(chimericRT.getJunctionPositions().contains(2100)); // already processed

        // assertEquals(1, chimericRT.getLocalChimericReads().size());
    }

}
