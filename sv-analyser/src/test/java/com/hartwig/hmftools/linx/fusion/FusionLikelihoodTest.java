package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.linx.gene.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.linx.fusion.FusionLikelihood.checkAddCombinedGenePhaseRegion;
import static com.hartwig.hmftools.linx.fusion.FusionLikelihood.calcOverlapBucketAreas;
import static com.hartwig.hmftools.linx.fusion.FusionLikelihood.generateGenePhaseRegions;
import static com.hartwig.hmftools.linx.fusion.GenePhaseRegion.PHASE_1;
import static com.hartwig.hmftools.linx.fusion.GenePhaseRegion.PHASE_2;
import static com.hartwig.hmftools.linx.fusion.GenePhaseRegion.PHASE_MAX;
import static com.hartwig.hmftools.linx.fusion.GenePhaseRegion.PHASE_NON_CODING;
import static com.hartwig.hmftools.linx.fusion.GenePhaseRegion.calcCombinedPhase;
import static com.hartwig.hmftools.linx.fusion.GenePhaseRegion.simpleToCombinedPhase;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import static junit.framework.TestCase.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptExonData;

import org.junit.Test;

public class FusionLikelihoodTest
{
    @Test
    public void testGeneRegions()
    {
        // 2 genes, each with transcripts with different phasings
        String geneId = "G001";
        byte strand = 1;
        EnsemblGeneData geneData = createEnsemblGeneData(geneId, "GEN1", "1", strand, 10, 100);

        List<TranscriptExonData> transExonDataList = Lists.newArrayList();

        // 3 coding exons, with coding region starting and ending half way through them

        String transName = "T001";
        int transId = 1;
        int transStart = 10;
        int transEnd = 80;
        Long codingStart = new Long(35);
        Long codingEnd = new Long(75);

        transExonDataList.add(new TranscriptExonData(geneId, transName, transId, true, strand, transStart, transEnd,
                10, 20, 1, -1, -1, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, transName, transId, true, strand, transStart, transEnd,
                30, 40, 2, -1, 1, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, transName, transId, true, strand, transStart, transEnd,
                50, 60, 3, 1, 2, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, transName, transId, true, strand, transStart, transEnd,
                70, 80, 4, 2, -1, codingStart, codingEnd, ""));

        transName = "T002";
        transId = 2;
        codingStart = new Long(15);
        codingEnd = new Long(55);


        transExonDataList.add(new TranscriptExonData(geneId, transName, transId, true, strand, transStart, transEnd,
                10, 20, 1, -1, 0, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, transName, transId, true, strand, transStart, transEnd,
                30, 40, 2, 0, 0, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, transName, transId, true, strand, transStart, transEnd,
                50, 60, 3, 0, -1, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, transName, transId, true, strand, transStart, transEnd,
                70, 80, 4, -1, -1, codingStart, codingEnd, ""));

        // and a non-coding transcript
        transName = "T003";
        transId = 3;
        codingStart = null;
        codingEnd = null;

        transExonDataList.add(new TranscriptExonData(geneId, transName, transId, true, strand, transStart, transEnd,
                10, 20, 1, -1, -1, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, transName, transId, true, strand, transStart, transEnd,
                70, 80, 2, -1, -1, codingStart, codingEnd, ""));

        List<GenePhaseRegion> phaseRegions = generateGenePhaseRegions(geneData, transExonDataList);
        assertEquals(5, phaseRegions.size());

        // test for the reverse strand
        geneId = "G002";
        strand = -1;
        geneData = createEnsemblGeneData(geneId, "GEN2", "1", strand, 10, 100);

        transExonDataList.clear();

        transName = "T004";
        transId = 4;
        transStart = 10;
        transEnd = 80;
        codingStart = new Long(35);
        codingEnd = new Long(75);

        transExonDataList.add(new TranscriptExonData(geneId, transName, transId, true, strand, transStart, transEnd,
                10, 20, 4, -1, -1, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, transName, transId, true, strand, transStart, transEnd,
                30, 40, 3, 1, -1, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, transName, transId, true, strand, transStart, transEnd,
                50, 60, 2, 2, 1, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, transName, transId, true, strand, transStart, transEnd,
                70, 80, 1, -1, 2, codingStart, codingEnd, ""));

        transName = "T005";
        transId = 5;
        codingStart = new Long(15);
        codingEnd = new Long(55);


        transExonDataList.add(new TranscriptExonData(geneId, transName, transId, true, strand, transStart, transEnd,
                10, 20, 4, 0, -1, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, transName, transId, true, strand, transStart, transEnd,
                30, 40, 3, 0, 0, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, transName, transId, true, strand, transStart, transEnd,
                50, 60, 2, -1, 0, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, transName, transId, true, strand, transStart, transEnd,
                70, 80, 1, -1, -1, codingStart, codingEnd, ""));

        // and a non-coding transcript
        transName = "T006";
        transId = 6;
        codingStart = null;
        codingEnd = null;

        transExonDataList.add(new TranscriptExonData(geneId, transName, transId, true, strand, transStart, transEnd,
                10, 20, 2, -1, -1, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, transName, transId, true, strand, transStart, transEnd,
                70, 80, 1, -1, -1, codingStart, codingEnd, ""));

        phaseRegions = generateGenePhaseRegions(geneData, transExonDataList);
        assertEquals(5, phaseRegions.size());
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

        assertTrue(hasPhaseRegion(regionsList, 50, 99, simpleToCombinedPhase(PHASE_2)));
        assertTrue(hasPhaseRegion(regionsList, 401, 499, simpleToCombinedPhase(PHASE_2)));
        assertTrue(hasPhaseRegion(regionsList, 601, 800, simpleToCombinedPhase(PHASE_2)));

        boolean[] phaseArray = new boolean[PHASE_MAX];
        phaseArray[PHASE_1] = true;
        phaseArray[PHASE_2] = true;
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

        GenePhaseRegion region = new GenePhaseRegion(geneId, 100, 200, PHASE_NON_CODING);
        assertTrue(region.hasPhaseOnly(PHASE_NON_CODING));

        region.addPhases(phaseArray);
        assertFalse(region.hasPhaseOnly(PHASE_NON_CODING));

        boolean[] phaseArray1 = {false, true, false, true, false};
        GenePhaseRegion region1 = new GenePhaseRegion(geneId, 100, 200, phaseArray1);

        boolean[] phaseArray2 = {true, false, true, false, true};
        GenePhaseRegion region2 = new GenePhaseRegion(geneId, 100, 200, phaseArray2);

        assertFalse(region1.hasAnyPhaseMatch(region2.getPhaseArray()));
        assertFalse(region2.hasAnyPhaseMatch(region1.getPhaseArray()));

        region2.addPhases(phaseArray1);
        region1.addPhases(phaseArray2);
        assertTrue(region1.hasAnyPhaseMatch(region2.getPhaseArray()));
        assertTrue(region2.hasAnyPhaseMatch(region1.getPhaseArray()));

        EnsemblGeneData gene = createEnsemblGeneData(geneId, "GEN2", "1", 1, 10, 100);
        GeneRangeData geneData = new GeneRangeData(gene);
        geneData.addPhaseRegions(Lists.newArrayList(region));
        assertFalse(geneData.hasCodingTranscripts());
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
    public void testProximateFusionCounts()
    {

        String geneId1 = "ESNG001";
        EnsemblGeneData gene1 = createEnsemblGeneData(geneId1, "GEN1", "1", 1, 10000, 12000);
        GeneRangeData lowerGene = new GeneRangeData(gene1);

        String geneId2 = "ESNG002";
        EnsemblGeneData gene2 = createEnsemblGeneData(geneId2, "GEN2", "1", 1, 10000, 12000);
        GeneRangeData upperGene = new GeneRangeData(gene1);

        List<Long> delLengths = Lists.newArrayList((long)50, (long)400);

        // first test 2 regions where the overlap is full within one of the del buckets (the second one)
        GenePhaseRegion lowerRegion = new GenePhaseRegion(geneId1, 100, 200, PHASE_1);
        GenePhaseRegion upperRegion = new GenePhaseRegion(geneId2, 300, 400, PHASE_1);

        Map<Integer, Long> bucketOverlapCounts = calcOverlapBucketAreas(delLengths, lowerGene, upperGene, lowerRegion, upperRegion, true);

        assertEquals(1, bucketOverlapCounts.size());
        long overlap = bucketOverlapCounts.get(0);
        assertEquals(10000, overlap);

        // with a max bucket length restriction
        delLengths.set(1, (long)150);

        bucketOverlapCounts = calcOverlapBucketAreas(delLengths, lowerGene, upperGene, lowerRegion, upperRegion, true);

        assertEquals(1, bucketOverlapCounts.size());
        overlap = bucketOverlapCounts.get(0);
        assertEquals(1275, overlap);

        // now with a min bucket length restriction
        delLengths.set(0, (long)250);
        delLengths.set(1, (long)400);

        bucketOverlapCounts = calcOverlapBucketAreas(delLengths, lowerGene, upperGene, lowerRegion, upperRegion, true);

        assertEquals(1, bucketOverlapCounts.size());
        overlap = bucketOverlapCounts.get(0);
        assertEquals(1275, overlap);

        // and with both
        delLengths.set(0, (long)225);
        delLengths.set(1, (long)275);

        bucketOverlapCounts = calcOverlapBucketAreas(delLengths, lowerGene, upperGene, lowerRegion, upperRegion, true);

        assertEquals(1, bucketOverlapCounts.size());
        overlap = bucketOverlapCounts.get(0);
        assertEquals(2525, overlap);
    }

    /*
    @Test
    public void testGeneRegionCounts()
    {
        String geneId  = "G001";
        byte strand = 1;
        EnsemblGeneData geneData = createEnsemblGeneData(geneId, "GEN1", "1", strand, 10, 100);

        List<TranscriptExonData> transExonDataList = Lists.newArrayList();

        // 3 coding exons, with coding region starting and ending half way through them

        int transStart = 10;
        int transEnd = 100;
        Long codingStart = new Long(35);
        Long codingEnd = new Long(75);

        transExonDataList.add(new TranscriptExonData(geneId, "T001", 1, true, strand, transStart, transEnd,
                10, 20, 1, -1, -1, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, "T001", 1, true, strand, transStart, transEnd,
                30, 40, 1, -1, 1, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, "T001", 1, true, strand, transStart, transEnd,
                50, 60, 1, 1, 2, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, "T001", 1, true, strand, transStart, transEnd,
                70, 80, 1, 2, -1, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, "T001", 1, true, strand, transStart, transEnd,
                90, 100, 1, -1, -1, codingStart, codingEnd, ""));

        int[] regionTotals5 = new int[GENE_PHASING_REGION_MAX];
        int[] regionTotals3 = new int[GENE_PHASING_REGION_MAX];
        setGenePhasingCounts(geneData, transExonDataList, regionTotals5, regionTotals3);

        assertEquals(regionTotals3[GENE_PHASING_REGION_5P_UTR], 26);
        //assertEquals(regionTotals[GENE_PHASING_REGION_CODING_0], 7);
        //assertEquals(regionTotals[GENE_PHASING_REGION_CODING_1], 17);
        //assertEquals(regionTotals[GENE_PHASING_REGION_CODING_2], 16);


        // now test with the reverse strand
        geneId  = "G002";
        strand = -1;
        EnsemblGeneData geneData2 = createEnsemblGeneData(geneId, "GEN2", "1", strand, 10, 100);

        transExonDataList.clear();

        // 3 coding exons, with coding region starting and ending half way through them

        codingStart = new Long(35);
        codingEnd = new Long(75);

        String transId = "T001";

        transExonDataList.add(new TranscriptExonData(geneId, transId, 1, true, strand, transStart, transEnd,
                10, 20, 1, -1, -1, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, transId, 1, true, strand, transStart, transEnd,
                30, 40, 1, -1, 1, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, transId, 1, true, strand, transStart, transEnd,
                50, 60, 1, 1, 2, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, transId, 1, true, strand, transStart, transEnd,
                70, 80, 1, 2, -1, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, transId, 1, true, strand, transStart, transEnd,
                90, 100, 1, -1, -1, codingStart, codingEnd, ""));

        setGenePhasingCounts(geneData, transExonDataList, regionTotals5, regionTotals3);

        // assertEquals(regionTotals[GENE_PHASING_REGION_5P_UTR], 26);
        //assertEquals(regionTotals[GENE_PHASING_REGION_CODING_0], 7);
        //assertEquals(regionTotals[GENE_PHASING_REGION_CODING_1], 17);
        //assertEquals(regionTotals[GENE_PHASING_REGION_CODING_2], 16);

    }
    */

}
