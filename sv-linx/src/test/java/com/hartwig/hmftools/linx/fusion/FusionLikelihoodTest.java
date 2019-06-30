package com.hartwig.hmftools.linx.fusion;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.linx.fusion_likelihood.CohortExpFusions.createPhaseRegionsFromTranscript;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseRegion.hasAnyPhaseMatch;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseRegion.hasNoOverlappingRegions;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseRegion.regionsPhaseMatched;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseType.PHASE_0;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseType.PHASE_1;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseType.PHASE_2;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseType.PHASE_5P_UTR;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseType.PHASE_MAX;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseType.PHASE_NON_CODING;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseType.typeAsInt;
import static com.hartwig.hmftools.linx.fusion_likelihood.LikelihoodCalc.calcOverlapBucketAreas;
import static com.hartwig.hmftools.linx.fusion_likelihood.PhaseRegionUtils.checkAddCombinedGenePhaseRegion;
import static com.hartwig.hmftools.linx.fusion_likelihood.PhaseRegionUtils.splitOverlappingPhaseRegion;
import static com.hartwig.hmftools.linx.gene.GeneTestUtils.addGeneData;
import static com.hartwig.hmftools.linx.gene.GeneTestUtils.addTransExonData;
import static com.hartwig.hmftools.linx.gene.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseRegion.calcCombinedPhase;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseRegion.simpleToCombinedPhase;
import static com.hartwig.hmftools.linx.gene.GeneTestUtils.createTransExons;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import static junit.framework.TestCase.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.linx.fusion_likelihood.CohortExpFusions;
import com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseRegion;
import com.hartwig.hmftools.linx.fusion_likelihood.GeneRangeData;
import com.hartwig.hmftools.linx.fusion_likelihood.RegionAllocator;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;

import org.junit.Test;

public class FusionLikelihoodTest
{
    @Test
    public void testGeneRegions()
    {
        SvGeneTranscriptCollection geneTransCache = new SvGeneTranscriptCollection();

        // 2 genes, each with transcripts with different phasings
        String geneId = "G001";
        byte strand = 1;
        EnsemblGeneData geneData = createEnsemblGeneData(geneId, "GEN1", "1", strand, 10, 100);

        GeneRangeData geneRangeData = new GeneRangeData(geneData);

        // List<TranscriptData> transcriptsList = Lists.newArrayList();
        List<TranscriptData> transDataList = Lists.newArrayList();

        // 3 coding exons, with coding region starting and ending half way through them

        int transId = 1;

        long[] exonStarts = new long[]{110, 130, 150, 170};
        int[] exonPhases = new int[]{-1, 1, 2, -1};
        TranscriptData transData = createTransExons(geneId, transId++, strand, exonStarts, exonPhases, 10);
        transDataList.add(transData);

        // converts to 110-129 -1, 130-149 1, 150-169 2

        // first test phase creation by itself
        List<GenePhaseRegion> phaseRegions = createPhaseRegionsFromTranscript(geneData, transData, 0);

        assertEquals(3, phaseRegions.size());
        assertTrue(hasPhaseRegion(phaseRegions, 110, 129, 10, 0));
        assertTrue(hasPhaseRegion(phaseRegions, 130, 149, 1000, 0));
        assertTrue(hasPhaseRegion(phaseRegions, 150, 169, 10000, 0));

        exonPhases = new int[]{0, 0, -1, -1};
        transData = createTransExons(geneId, transId++, strand, exonStarts, exonPhases, 10);
        transDataList.add(transData);

        // converts to 110-129 0, 130-149 0

        phaseRegions = createPhaseRegionsFromTranscript(geneData, transData, 0);

        assertEquals(2, phaseRegions.size());
        assertTrue(hasPhaseRegion(phaseRegions, 110, 129, 100, 0));
        assertTrue(hasPhaseRegion(phaseRegions, 130, 149, 100, 0));

        // and a non-coding transcript
        exonStarts = new long[]{110, 170};
        exonPhases = new int[]{-1, -1};
        transData = createTransExons(geneId, transId++, strand, exonStarts, exonPhases, 10);
        transDataList.add(transData);

        // converts to 110-179 -1

        phaseRegions = createPhaseRegionsFromTranscript(geneData, transData, 0);

        assertEquals(1, phaseRegions.size());
        assertTrue(hasPhaseRegion(phaseRegions, 110, 180, 1, 0));

        // now test assigning a set of transcripts to a gene
        CohortExpFusions likelihoodCalc = new CohortExpFusions();

        // for testing same-gene fusions
        List<Long> delLengths = Lists.newArrayList((long)1, (long)1000);
        likelihoodCalc.initialise(delLengths, 0);

        likelihoodCalc.generatePhaseRegions(geneRangeData, transDataList, geneTransCache);
        phaseRegions = geneRangeData.getPhaseRegions();
        assertEquals(4, phaseRegions.size());

        assertTrue(hasPhaseRegion(phaseRegions, 110, 129, 111, 0));
        assertTrue(hasPhaseRegion(phaseRegions, 130, 149, 1101, 0));
        assertTrue(hasPhaseRegion(phaseRegions, 150, 169, 10001, 0));
        assertTrue(hasPhaseRegion(phaseRegions, 170, 180, 1, 0));

        Map<Integer,Long> transSAMap = geneTransCache.getTransSpliceAcceptorPosDataMap();

        transDataList.stream().forEach(x -> transSAMap.put(x.TransId, (long)60));

        // test again but with a preceding gene region distance
        likelihoodCalc.generatePhaseRegions(geneRangeData, transDataList, geneTransCache);
        phaseRegions = geneRangeData.getPhaseRegions();
        assertEquals(5, phaseRegions.size());

        assertTrue(hasPhaseRegion(phaseRegions, 60, 109, 110, 110));

        // test for the reverse strand
        geneId = "G002";
        strand = -1;
        geneData = createEnsemblGeneData(geneId, "GEN2", "1", strand, 10, 100);

        geneRangeData = new GeneRangeData(geneData);

        transDataList.clear();

        exonStarts = new long[]{10, 30, 50, 70};
        exonPhases = new int[]{-1, -1, 1, 2};
        transData = createTransExons(geneId, transId++, strand, exonStarts, exonPhases, 10);
        transDataList.add(transData);

        // converts to 41-60 1, 61-80 2

        exonStarts = new long[]{10, 30, 50, 70};
        exonPhases = new int[]{-1, 0, 0, -1};
        transData = createTransExons(geneId, transId++, strand, exonStarts, exonPhases, 10);
        transDataList.add(transData);

        // converts to 21-40 0, 41-60 0, 61-80 -1

        // and a non-coding transcript
        exonStarts = new long[]{10, 70};
        exonPhases = new int[]{-1, -1};

        // converts to 20-80 -1 NC

        transData = createTransExons(geneId, transId++, strand, exonStarts, exonPhases, 10);
        transDataList.add(transData);

        transSAMap.clear();

        likelihoodCalc.generatePhaseRegions(geneRangeData, transDataList, geneTransCache);
        phaseRegions = geneRangeData.getPhaseRegions();
        assertEquals(4, phaseRegions.size());

        assertTrue(hasPhaseRegion(phaseRegions, 10, 20, 1, 0));
        assertTrue(hasPhaseRegion(phaseRegions, 21, 40, 101, 0));
        assertTrue(hasPhaseRegion(phaseRegions, 41, 60, 1101, 0));
        assertTrue(hasPhaseRegion(phaseRegions, 61, 80, 10011, 0));

        transDataList.stream().forEach(x -> transSAMap.put(x.TransId, (long)100000));

        // test again but with a preceding gene region distance - will be capped at 10K
        likelihoodCalc.generatePhaseRegions(geneRangeData, transDataList, geneTransCache);
        phaseRegions = geneRangeData.getPhaseRegions();
        assertEquals(5, phaseRegions.size());

        assertTrue(hasPhaseRegion(phaseRegions, 81, 10080, 10010, 10010));

        // now test with 3 transcripts each with the same exons but kept separate for same-gene tests
        transDataList.clear();

        geneId = "G003";
        strand = 1;
        geneData = createEnsemblGeneData(geneId, geneId, "1", strand, 10, 100);

        geneRangeData = new GeneRangeData(geneData);

        exonStarts = new long[]{10, 30, 50, 70, 90};
        exonPhases = new int[]{-1, 1, 1, 1, -1};
        transData = createTransExons(geneId, transId++, strand, exonStarts, exonPhases, 10);
        transDataList.add(transData);

        exonStarts = new long[]{10, 30, 50, 70, 90};
        exonPhases = new int[]{-1, 1, 1, 1, -1};
        transData = createTransExons(geneId, transId++, strand, exonStarts, exonPhases, 10);
        transDataList.add(transData);

        transSAMap.clear();

        transData = createTransExons(geneId, transId++, strand, exonStarts, exonPhases, 10);
        transDataList.add(transData);

        likelihoodCalc.generatePhaseRegions(geneRangeData, transDataList, geneTransCache);
        // likelihoodCalc.generateSameGeneCounts(geneRangeData);

        geneRangeData.clearOverlapCounts();

        /*
        likelihoodCalc.generatePhaseRegions(geneRangeData, transExonDataList, geneTransCache);
        assertEquals(1, geneRangeData.getDelFusionBaseCounts().size());

        Long overlapCount = new Long(19*19 +19*19 + 19*19);
        assertEquals(overlapCount, geneRangeData.getDelFusionBaseCounts().get(0));

        // test again but with no phase matches within the same transcript
        transExonDataList.clear();
        geneRangeData.clearOverlapCounts();
        exonStarts = new long[]{10, 30, 50, 70, 90};
        exonPhases = new int[]{-1, 0, 1, 2, -1};
        createTransExons(transExonDataList, geneId, transId++, strand, exonStarts, exonPhases, 10);

        exonStarts = new long[]{10, 30, 50, 70, 90};
        exonPhases = new int[]{-1, 2, 1, 0, -1};
        createTransExons(transExonDataList, geneId, transId++, strand, exonStarts, exonPhases, 10);

        likelihoodCalc.generatePhaseRegions(geneRangeData, transExonDataList, geneTransCache);
        assertTrue(geneRangeData.getDelFusionBaseCounts().isEmpty());
        */
    }

    private static boolean hasPhaseRegion(List<GenePhaseRegion> regionsList, long start, long end, int combinedPhase, int combinedPreGene)
    {
        for(GenePhaseRegion region : regionsList)
        {
            if(region.start() == start && region.end() == end
            && region.getCombinedPhase() == combinedPhase && region.getCombinedPreGeneStatus() == combinedPreGene)
            {
                return true;
            }
        }

        return false;
    }


    @Test
    public void testPhaseRegionMerging()
    {
        final String geneId = "001";

        List<GenePhaseRegion> regionsList = Lists.newArrayList();

        // test enclosure and overlap for matching phases
        checkAddCombinedGenePhaseRegion(new GenePhaseRegion(geneId, 100, 400, PHASE_1), regionsList);
        checkAddCombinedGenePhaseRegion(new GenePhaseRegion(geneId, 200, 300, PHASE_1), regionsList);
        assertEquals(1, regionsList.size());

        assertTrue(hasPhaseRegion(regionsList, 100, 400, simpleToCombinedPhase(PHASE_1)));

        // add a non-overlapping phase
        checkAddCombinedGenePhaseRegion(new GenePhaseRegion(geneId, 500, 600, PHASE_1), regionsList);
        assertEquals(2, regionsList.size());

        // and now one which overlaps them both with a different phase
        checkAddCombinedGenePhaseRegion(new GenePhaseRegion(geneId, 50, 800, PHASE_2), regionsList);
        assertEquals(5, regionsList.size());
        assertTrue(hasNoOverlappingRegions(regionsList));

        assertTrue(hasPhaseRegion(regionsList, 50, 99, simpleToCombinedPhase(PHASE_2)));
        assertTrue(hasPhaseRegion(regionsList, 401, 499, simpleToCombinedPhase(PHASE_2)));
        assertTrue(hasPhaseRegion(regionsList, 601, 800, simpleToCombinedPhase(PHASE_2)));

        boolean[] phaseArray = new boolean[PHASE_MAX];
        boolean[] preGeneStatus = new boolean[PHASE_MAX];
        phaseArray[typeAsInt(PHASE_1)] = true;
        phaseArray[typeAsInt(PHASE_2)] = true;
        int combinedPhase = calcCombinedPhase(phaseArray);
        assertTrue(hasPhaseRegion(regionsList, 100, 400, combinedPhase));
        assertTrue(hasPhaseRegion(regionsList, 500, 600, combinedPhase));

        // test overlapping phases
        regionsList.clear();

        checkAddCombinedGenePhaseRegion(new GenePhaseRegion(geneId, 100, 300, PHASE_2), regionsList);
        checkAddCombinedGenePhaseRegion(new GenePhaseRegion(geneId, 200, 400, PHASE_1), regionsList);

        assertEquals(3, regionsList.size());
        assertTrue(hasPhaseRegion(regionsList, 100, 199, simpleToCombinedPhase(PHASE_2)));
        assertTrue(hasPhaseRegion(regionsList, 200, 300, combinedPhase));
        assertTrue(hasPhaseRegion(regionsList, 301, 400, simpleToCombinedPhase(PHASE_1)));
        assertTrue(hasNoOverlappingRegions(regionsList));

        GenePhaseRegion region = new GenePhaseRegion(geneId, 100, 200, PHASE_NON_CODING);
        assertTrue(region.hasPhaseOnly(PHASE_NON_CODING));

        region.addPhases(phaseArray, preGeneStatus);
        assertFalse(region.hasPhaseOnly(PHASE_NON_CODING));

        boolean[] phaseArray1 = {false, true, false, true, false};
        GenePhaseRegion region1 = new GenePhaseRegion(geneId, 100, 200, phaseArray1, preGeneStatus);

        boolean[] phaseArray2 = {true, false, true, false, true};
        GenePhaseRegion region2 = new GenePhaseRegion(geneId, 100, 200, phaseArray2, preGeneStatus);

        assertFalse(region1.hasAnyPhaseMatch(region2.getPhaseArray()));
        assertFalse(region2.hasAnyPhaseMatch(region1.getPhaseArray()));

        region2.addPhases(phaseArray1, preGeneStatus);
        region1.addPhases(phaseArray2, preGeneStatus);
        assertTrue(region1.hasAnyPhaseMatch(region2.getPhaseArray()));
        assertTrue(region2.hasAnyPhaseMatch(region1.getPhaseArray()));

        EnsemblGeneData gene = createEnsemblGeneData(geneId, "GEN2", "1", 1, 10, 100);
        GeneRangeData geneData = new GeneRangeData(gene);
        geneData.setPhaseRegions(Lists.newArrayList(region));
        assertFalse(geneData.hasCodingTranscripts());

        // test not expanding an existing region to overlap others just because it matches a new one exactly
        regionsList.clear();

        checkAddCombinedGenePhaseRegion(new GenePhaseRegion(geneId, 100, 200, PHASE_2), regionsList);
        checkAddCombinedGenePhaseRegion(new GenePhaseRegion(geneId, 201, 300, PHASE_1), regionsList);
        checkAddCombinedGenePhaseRegion(new GenePhaseRegion(geneId, 150, 250, PHASE_2), regionsList);

        assertEquals(3, regionsList.size());
        assertTrue(hasPhaseRegion(regionsList, 100, 200, 10000, 0));
        assertTrue(hasPhaseRegion(regionsList, 201, 250, 11000, 0));
        assertTrue(hasPhaseRegion(regionsList, 251, 300, 1000, 0));
    }

    @Test
    public void testPhaseSplitting()
    {
        GenePhaseRegion region1 = new GenePhaseRegion("G1", 100, 300, PHASE_0);
        GenePhaseRegion region2 = new GenePhaseRegion("G2", 200, 400, PHASE_0);

        List<GenePhaseRegion> regions1 = Lists.newArrayList(region1);
        List<GenePhaseRegion> regions2 = Lists.newArrayList(region2);

        splitOverlappingPhaseRegion(region1, 0, regions1, region2, 0, regions2);

        assertEquals(100, region1.start());
        assertEquals(250, region1.end());
        assertEquals(1, regions1.size());
        assertEquals(251, region2.start());
        assertEquals(400, region2.end());
        assertEquals(1, regions2.size());

        region1 = new GenePhaseRegion("G1", 200, 400, PHASE_0);
        region2 = new GenePhaseRegion("G2", 100, 300, PHASE_0);

        splitOverlappingPhaseRegion(region1, 0, regions1, region2, 0, regions2);

        assertEquals(251, region1.start());
        assertEquals(400, region1.end());
        assertEquals(1, regions1.size());
        assertEquals(100, region2.start());
        assertEquals(250, region2.end());
        assertEquals(1, regions2.size());

        // region 1 protein-coding, the other not
        region1 = new GenePhaseRegion("G1", 100, 300, PHASE_0);
        region2 = new GenePhaseRegion("G2", 200, 400, PHASE_0);
        region1.setProteinCoding(true);

        splitOverlappingPhaseRegion(region1, 0, regions1, region2, 0, regions2);

        assertEquals(100, region1.start());
        assertEquals(300, region1.end());
        assertEquals(1, regions1.size());
        assertEquals(301, region2.start());
        assertEquals(400, region2.end());
        assertEquals(1, regions2.size());

        // one enclosing the other
        region1 = new GenePhaseRegion("G1", 100, 400, PHASE_0);
        region2 = new GenePhaseRegion("G2", 200, 300, PHASE_0);

        splitOverlappingPhaseRegion(region1, 0, regions1, region2, 0, regions2);

        assertEquals(100, region1.start());
        assertEquals(250, region1.end());
        assertEquals(2, regions1.size());
        assertEquals(301, regions1.get(1).start());
        assertEquals(400, regions1.get(1).end());

        assertEquals(251, region2.start());
        assertEquals(300, region2.end());
        assertEquals(1, regions2.size());

        region1 = new GenePhaseRegion("G1", 100, 400, PHASE_0);
        region2 = new GenePhaseRegion("G2", 200, 300, PHASE_0);
        region1.setProteinCoding(true);
        regions1 = Lists.newArrayList(region1);
        regions2 = Lists.newArrayList(region2);

        splitOverlappingPhaseRegion(region1, 0, regions1, region2, 0, regions2);

        assertEquals(100, region1.start());
        assertEquals(400, region1.end());
        assertTrue( regions2.isEmpty());

        region1 = new GenePhaseRegion("G1", 100, 400, PHASE_0);
        region2 = new GenePhaseRegion("G2", 200, 300, PHASE_0);
        region2.setProteinCoding(true);
        regions1 = Lists.newArrayList(region1);
        regions2 = Lists.newArrayList(region2);

        splitOverlappingPhaseRegion(region1, 0, regions1, region2, 0, regions2);
        assertEquals(100, region1.start());
        assertEquals(199, region1.end());
        assertEquals(2, regions1.size());
        assertEquals(301, regions1.get(1).start());
        assertEquals(400, regions1.get(1).end());

        assertEquals(200, region2.start());
        assertEquals(300, region2.end());
        assertEquals(1, regions2.size());
    }

    @Test
    public void testPhaseMatching()
    {
        final String geneId = "001";

        GenePhaseRegion region1 = new GenePhaseRegion(geneId, 100, 200, PHASE_0);
        region1.setPreGene(true, region1.Phase);

        GenePhaseRegion region2 = new GenePhaseRegion(geneId, 100, 200, PHASE_0);

        assertFalse(hasAnyPhaseMatch(region1, region2, false));
        assertFalse(hasAnyPhaseMatch(region1, region2, true));
        assertTrue(hasAnyPhaseMatch(region2, region1, true));

        assertFalse(hasAnyPhaseMatch(region1, region2, false));
        assertTrue(regionsPhaseMatched(region2, region1));
        assertFalse(regionsPhaseMatched(region1, region2));

        // test with a non-coding 5' region to downstream UTR
        region1 = new GenePhaseRegion(geneId, 100, 200, PHASE_NON_CODING);

        region2 = new GenePhaseRegion(geneId, 100, 200, PHASE_5P_UTR);
        assertFalse(regionsPhaseMatched(region2, region1));
        assertTrue(regionsPhaseMatched(region1, region2));

        // don't allow standard non-coding to non-coding
        region2 = new GenePhaseRegion(geneId, 100, 200, PHASE_NON_CODING);
        assertFalse(hasAnyPhaseMatch(region1, region2, false));

        // test non-coding to downstream 5'UTR
        region2 = new GenePhaseRegion(geneId, 100, 200, PHASE_5P_UTR);
        region2.setPreGene(true, region2.Phase);

        assertFalse(hasAnyPhaseMatch(region1, region2, false));
        assertFalse(hasAnyPhaseMatch(region2, region1, false));
        assertFalse(hasAnyPhaseMatch(region1, region2, true));
        assertFalse(regionsPhaseMatched(region2, region1));
        assertTrue(regionsPhaseMatched(region1, region2));
    }

    private static boolean hasPhaseRegion(List<GenePhaseRegion> regionsList, long start, long end, int combinedPhase)
    {
        for(GenePhaseRegion region : regionsList)
        {
            if(region.start() == start && region.end() == end && region.getCombinedPhase() == combinedPhase)
                return true;
        }

        return false;
    }

    @Test
    public void testRegionAllocator()
    {
        RegionAllocator regionAllocator = new RegionAllocator(100);

        assertEquals(0, regionAllocator.baseToIndex(0));
        assertEquals(0, regionAllocator.baseToIndex(99));
        assertEquals(1, regionAllocator.baseToIndex(100));

        int blockArea = regionAllocator.blockSize() * regionAllocator.blockSize();

        // first test without any restricting bucket lengths
        long overlap = regionAllocator.allocateBases(10, 45, 55, 99);
        assertEquals(blockArea, overlap);

        // cannot reallocate
        overlap = regionAllocator.allocateBases(10, 45, 55, 99);
        assertEquals(0, overlap);

        regionAllocator.reset();

        overlap = regionAllocator.allocateBases(0, 499, 800, 999);
        assertEquals(blockArea * 5 * 2, overlap);

        overlap = regionAllocator.allocateBases(300, 599, 600, 999);
        assertEquals(blockArea * (3 * 2 + 2), overlap);

        regionAllocator.reset();

        // now test with length min and max
        int minBucketLen = 100;
        int maxBucketLen = 300;
        overlap = regionAllocator.allocateBases(100, 599, 100, 599, minBucketLen, maxBucketLen, true);
        assertEquals(blockArea * (3 + 3 + 2 + 1), overlap);

    }

    /*
    @Test
    public void testSameGeneCounts()
    {
        EnsemblGeneData gene = createEnsemblGeneData("ESNG001", "GEN1", "1", 1, 10000, 12000);
        GeneRangeData geneData = new GeneRangeData(gene);

        List<Long> delLengths = Lists.newArrayList((long)1, (long)1000);

        List<GenePhaseRegion> transcriptRegions = Lists.newArrayList();

        transcriptRegions.add(new GenePhaseRegion(gene.GeneId, 100, 199, PHASE_0));
        transcriptRegions.add(new GenePhaseRegion(gene.GeneId, 300, 399, PHASE_0));

        for(int i = 0; i < transcriptRegions.size(); ++i)
        {
            transcriptRegions.get(i).setTransId(0);
        }

        // second set of transcripts matches exactly and so won't be counted twice
        transcriptRegions.add(new GenePhaseRegion(gene.GeneId, 100, 199, PHASE_0));
        transcriptRegions.add(new GenePhaseRegion(gene.GeneId, 300, 399, PHASE_0));

        for(int i = 2; i < transcriptRegions.size(); ++i)
        {
            transcriptRegions.get(i).setTransId(1);
        }

        // a third set which overlaps in part with the first set
        transcriptRegions.add(new GenePhaseRegion(gene.GeneId, 100, 199, PHASE_0));
        transcriptRegions.add(new GenePhaseRegion(gene.GeneId, 300, 399, PHASE_0));
        transcriptRegions.add(new GenePhaseRegion(gene.GeneId, 500, 599, PHASE_0));

        for(int i = 4; i < transcriptRegions.size(); ++i)
        {
            transcriptRegions.get(i).setTransId(2);
        }

        geneData.setTranscriptPhaseRegions(transcriptRegions);

        CohortExpFusions likelihoodCalc = new CohortExpFusions();
        likelihoodCalc.initialise(delLengths, 0);
        likelihoodCalc.generateSameGeneCounts(geneData);

        Long overlapCount = new Long((100 * 100) * 3);

        assertEquals(1, geneData.getDelFusionBaseCounts().size());
        assertEquals(overlapCount, geneData.getDelFusionBaseCounts().get(0));
    }
    */

    @Test
    public void testProximateFusionCounts()
    {
        EnsemblGeneData gene1 = createEnsemblGeneData("ESNG001", "GEN1", "1", 1, 10000, 12000);
        GeneRangeData lowerGene = new GeneRangeData(gene1);

        EnsemblGeneData gene2 = createEnsemblGeneData("ESNG002", "GEN2", "1", 1, 10000, 12000);
        GeneRangeData upperGene = new GeneRangeData(gene2);

        List<Long> delLengths = Lists.newArrayList((long)50, (long)400);

        // first test 2 regions where the overlap is full within one of the del buckets (the second one)
        GenePhaseRegion lowerRegion = new GenePhaseRegion(gene1.GeneId, 100, 200, PHASE_1);
        GenePhaseRegion upperRegion = new GenePhaseRegion(gene2.GeneId, 300, 400, PHASE_1);

        Map<Integer, Long> bucketOverlapCounts = calcOverlapBucketAreas(delLengths, null, lowerGene, upperGene, lowerRegion, upperRegion, true);

        assertEquals(1, bucketOverlapCounts.size());
        long overlap = bucketOverlapCounts.get(0);
        assertEquals(10000, overlap);

        // with a max bucket length restriction
        delLengths.set(1, (long)150);

        bucketOverlapCounts = calcOverlapBucketAreas(delLengths, null, lowerGene, upperGene, lowerRegion, upperRegion, true);

        assertEquals(1, bucketOverlapCounts.size());
        overlap = bucketOverlapCounts.get(0);
        assertEquals(1275, overlap);

        // now with a min bucket length restriction
        delLengths.set(0, (long)250);
        delLengths.set(1, (long)400);

        bucketOverlapCounts = calcOverlapBucketAreas(delLengths, null, lowerGene, upperGene, lowerRegion, upperRegion, true);

        assertEquals(1, bucketOverlapCounts.size());
        overlap = bucketOverlapCounts.get(0);
        assertEquals(1275, overlap);

        // and with both
        delLengths.set(0, (long)225);
        delLengths.set(1, (long)275);

        bucketOverlapCounts = calcOverlapBucketAreas(delLengths, null, lowerGene, upperGene, lowerRegion, upperRegion, true);

        assertEquals(1, bucketOverlapCounts.size());
        overlap = bucketOverlapCounts.get(0);
        assertEquals(2525, overlap);
    }

    /*
    @Test
    public void testOverlappingProximateRegions()
    {
        List<Long> delDupLengths = Lists.newArrayList((long)10, (long)10000);

        CohortExpFusions likelihoodCalc = new CohortExpFusions();
        likelihoodCalc.initialise(delDupLengths, 0);

        GeneRangeData lowerGene = new GeneRangeData(createEnsemblGeneData("G1", "G1", "1", 1, 100, 1000));
        GeneRangeData upperGene = new GeneRangeData(createEnsemblGeneData("G2", "G2", "1", 1, 100, 1000));

        GenePhaseRegion lowerRegion = new GenePhaseRegion("G1", 100, 300, PHASE_0);
        GenePhaseRegion upperRegion = new GenePhaseRegion("G2", 200, 400, PHASE_0);
        likelihoodCalc.testOverlappingProximateRegions(lowerGene, upperGene, lowerRegion, upperRegion);

        assertEquals(1, lowerGene.getDelFusionBaseCounts().size());
        assertEquals(1, upperGene.getDelFusionBaseCounts().size());

        long overlap = (long)(200 * 200);

        assertTrue(approxEqual(overlap,  lowerGene.getDelFusionBaseCounts().get(0), 0.1));

    }
    */

    /*
    @Test
    public void testIntegratedCounts()
    {
        // test multiple chromosomes together for the various types of counts

        SvGeneTranscriptCollection geneTransCache = new SvGeneTranscriptCollection();

        int geneIndex = 1;
        String geneId, geneName;

        int transId = 1;

        List<EnsemblGeneData> geneList = Lists.newArrayList();

        // for 3 chromosomes add 3 genes on one strand, 3 on the other at varying distances
        // for calculate simplicity make each intron distance 1000 bases

        long[] geneExons1 = {1000, 2000, 3000, 4000, 5000};
        long[] geneExons2 = {10000, 12000, 13000};
        long[] geneExons3 = {200000, 202000, 203000};

        for(int i = 1; i <= 3; ++i)
        {
            String chromosome = String.valueOf(i);

            for(int j = 0; j <= 1; ++j)
            {
                byte strand = (j == 0) ? (byte)1 : (byte)-1;

                geneId = geneName = String.format("ESNG%04d", geneIndex++);
                TranscriptData transData = createTransExons(geneId, transId++, strand, geneExons1, new int[] {-1, 0, 1, 0, -1}, 100);
                addTransExonData(geneTransCache, geneId, Lists.newArrayList(transData));
                geneList.add(createEnsemblGeneData(geneId, geneName, chromosome, strand, transData.TransStart, transData.TransEnd));

                geneId = geneName = String.format("ESNG%04d", geneIndex++);
                transData = createTransExons(geneId, transId++, strand, geneExons2, new int[] {-1, 0, -1}, 100);
                addTransExonData(geneTransCache, geneId, Lists.newArrayList(transData));
                geneList.add(createEnsemblGeneData(geneId, geneName, chromosome, strand, transData.TransStart, transData.TransEnd));

                // another transcript outside the DEL and DUP bucket length ranges
                geneId = geneName = String.format("ESNG%04d", geneIndex++);
                transData = createTransExons(geneId, transId++, strand, geneExons3, new int[] {-1, 0, -1}, 100);
                addTransExonData(geneTransCache, geneId, Lists.newArrayList(transData));
                geneList.add(createEnsemblGeneData(geneId, geneName, chromosome, strand, transData.TransStart, transData.TransEnd));
            }

            addGeneData(geneTransCache, chromosome, geneList);
            geneList = Lists.newArrayList();
        }

        // DELs can link sae gene and both closer genes, DUPs on the closer genes
        List<Long> delDupLengths = Lists.newArrayList((long)50, (long)5000, (long)50000);
        int shortInv = 20000;


        CohortExpFusions likelihoodCalc = new CohortExpFusions();
        likelihoodCalc.initialise(delDupLengths, shortInv);
        likelihoodCalc.generateExpectedFusions(geneTransCache, Lists.newArrayList(), Lists.newArrayList());

        // expected overlap counts
        int intronLength = 1000;
        long shortDelOverlap = intronLength * intronLength;
        long medDelOverlap = 2 * intronLength * intronLength;

        final Map<String, List<GeneRangeData>> chrGeneDataMap = likelihoodCalc.getChrGeneRangeDataMap();
        assertEquals(3, chrGeneDataMap.size());

        List<GeneRangeData> geneRangeList = chrGeneDataMap.get("1");
        assertEquals(6, geneRangeList.size());

        GeneRangeData geneData = geneRangeList.get(0);
        assertEquals(4, geneData.getPhaseRegions().size());
        assertEquals(4, geneData.getPhaseRegions().size());

        assertEquals(2, geneData.getDelFusionBaseCounts().size());

        assertTrue(approxEqual(shortDelOverlap, geneData.getDelFusionBaseCounts().get(0), 0.01));
    }
    */

    private boolean approxEqual(long bases1, long bases2, double allowPerc)
    {
        if((bases1 == 0) != (bases2 == 0))
            return false;

        double diff1 = (abs(bases1 - bases2) / (double)bases1);
        double diff2 = (abs(bases1 - bases2) / (double)bases2);

        return diff1 <= allowPerc && diff2 <= allowPerc;
    }

}
