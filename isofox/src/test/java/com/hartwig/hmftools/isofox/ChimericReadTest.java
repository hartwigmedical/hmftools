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
import com.hartwig.hmftools.isofox.fusion.FusionFragment;

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

        // a DEL more than 500K from the gene boundary is considered chimeric
        chimericRT.clear();

        // DEL linking 2 genes at known splice sites is considered chimeric
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

}
