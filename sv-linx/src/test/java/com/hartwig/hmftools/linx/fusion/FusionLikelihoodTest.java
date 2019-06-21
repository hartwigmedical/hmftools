package com.hartwig.hmftools.linx.fusion;

import static java.lang.Math.abs;

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
import static com.hartwig.hmftools.linx.fusion_likelihood.PhaseRegionUtils.validateSimpleVsCombinedPhaseRegions;
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
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptExonData;
import com.hartwig.hmftools.linx.fusion_likelihood.CohortExpFusions;
import com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseRegion;
import com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseType;
import com.hartwig.hmftools.linx.fusion_likelihood.GeneRangeData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;

import org.junit.Test;

public class FusionLikelihoodTest
{
    @Test
    public void testGeneRegions()
    {
        CohortExpFusions likelihoodCalc = new CohortExpFusions();

        // for testing same-gene fusions
        List<Long> delLengths = Lists.newArrayList((long)1, (long)1000);
        likelihoodCalc.initialise(delLengths, 0, false);

        // 2 genes, each with transcripts with different phasings
        String geneId = "G001";
        byte strand = 1;
        EnsemblGeneData geneData = createEnsemblGeneData(geneId, "GEN1", "1", strand, 10, 100);

        GeneRangeData geneRangeData = new GeneRangeData(geneData);

        List<TranscriptExonData> transExonDataList = Lists.newArrayList();

        // 3 coding exons, with coding region starting and ending half way through them

        int transId = 1;

        long[] exonStarts = new long[]{110, 130, 150, 170};
        int[] exonPhases = new int[]{-1, 1, 2, -1};
        createTransExons(transExonDataList, geneId, transId++, strand, exonStarts, exonPhases, 10);

        exonPhases = new int[]{0,0, -1, -1};
        createTransExons(transExonDataList, geneId, transId++, strand, exonStarts, exonPhases, 10);

        // and a non-coding transcript
        exonStarts = new long[]{110, 170};
        exonPhases = new int[]{-1, -1};
        createTransExons(transExonDataList, geneId, transId++, strand, exonStarts, exonPhases, 10);

        likelihoodCalc.generateGenePhaseRegions(geneRangeData, transExonDataList, 0);
        List<GenePhaseRegion> phaseRegions = geneRangeData.getPhaseRegions();
        assertEquals(5, phaseRegions.size());

        assertTrue(hasPhaseRegion(phaseRegions, 110, 129, PHASE_5P_UTR, false));
        assertTrue(hasPhaseRegion(phaseRegions, 110, 149, PHASE_0, false));
        assertTrue(hasPhaseRegion(phaseRegions, 130, 149, PHASE_1, false));
        assertTrue(hasPhaseRegion(phaseRegions, 150, 169, PHASE_2, false));
        assertTrue(hasPhaseRegion(phaseRegions, 110, 180, PHASE_NON_CODING, false));

        // test again but with a preceding gene region distance
        likelihoodCalc.generateGenePhaseRegions(geneRangeData, transExonDataList, 60);
        phaseRegions = geneRangeData.getPhaseRegions();
        assertEquals(7, phaseRegions.size());

        assertTrue(hasPhaseRegion(phaseRegions, 60, 109, PHASE_5P_UTR, true));
        assertTrue(hasPhaseRegion(phaseRegions, 60, 109, PHASE_0, true));

        // test for the reverse strand
        geneId = "G002";
        strand = -1;
        geneData = createEnsemblGeneData(geneId, "GEN2", "1", strand, 10, 100);

        geneRangeData = new GeneRangeData(geneData);

        transExonDataList.clear();

        exonStarts = new long[]{10, 30, 50, 70};
        exonPhases = new int[]{-1, -1, 1, 2};
        createTransExons(transExonDataList, geneId, transId++, strand, exonStarts, exonPhases, 10);

        exonStarts = new long[]{10, 30, 50, 70};
        exonPhases = new int[]{-1, 0, 0, -1};
        createTransExons(transExonDataList, geneId, transId++, strand, exonStarts, exonPhases, 10);

        // and a non-coding transcript
        exonStarts = new long[]{10, 70};
        exonPhases = new int[]{-1, -1};

        createTransExons(transExonDataList, geneId, transId++, strand, exonStarts, exonPhases, 10);

        likelihoodCalc.generateGenePhaseRegions(geneRangeData, transExonDataList, 0);
        phaseRegions = geneRangeData.getPhaseRegions();
        assertEquals(5, phaseRegions.size());

        assertTrue(hasPhaseRegion(phaseRegions, 61, 80, PHASE_5P_UTR, false));
        assertTrue(hasPhaseRegion(phaseRegions, 21, 60, PHASE_0, false));
        assertTrue(hasPhaseRegion(phaseRegions, 41, 60, PHASE_1, false));
        assertTrue(hasPhaseRegion(phaseRegions, 61, 80, PHASE_2, false));
        assertTrue(hasPhaseRegion(phaseRegions, 10, 80, PHASE_NON_CODING, false));

        // test again but with a preceding gene region distance - will be capped at 10K
        likelihoodCalc.generateGenePhaseRegions(geneRangeData, transExonDataList, 100000);
        phaseRegions = geneRangeData.getPhaseRegions();
        assertEquals(7, phaseRegions.size());

        assertTrue(hasPhaseRegion(phaseRegions, 81, 10080, PHASE_5P_UTR, true));
        assertTrue(hasPhaseRegion(phaseRegions, 81, 10080, PHASE_2, true));

        // now test with 3 transcripts each with the same exons but kept seperate for same-gene tests
        transExonDataList.clear();

        geneId = "G003";
        strand = 1;
        geneData = createEnsemblGeneData(geneId, geneId, "1", strand, 10, 100);

        geneRangeData = new GeneRangeData(geneData);

        exonStarts = new long[]{10, 30, 50, 70, 90};
        exonPhases = new int[]{-1, 1, 1, 1, -1};
        createTransExons(transExonDataList, geneId, transId++, strand, exonStarts, exonPhases, 10);

        exonStarts = new long[]{10, 30, 50, 70, 90};
        exonPhases = new int[]{-1, 1, 1, 1, -1};
        createTransExons(transExonDataList, geneId, transId++, strand, exonStarts, exonPhases, 10);

        geneRangeData.clearOverlapCounts();
        likelihoodCalc.generateGenePhaseRegions(geneRangeData, transExonDataList, 0);
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

        likelihoodCalc.generateGenePhaseRegions(geneRangeData, transExonDataList, 0);
        assertTrue(geneRangeData.getDelFusionBaseCounts().isEmpty());
    }

    private static boolean hasPhaseRegion(List<GenePhaseRegion> regionsList, long start, long end, GenePhaseType phase, boolean isPreGene)
    {
        for(GenePhaseRegion region : regionsList)
        {
            if(region.start() == start && region.end() == end && region.Phase == phase && region.isAnyPreGene() == isPreGene)
                return true;
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
    }

    @Test
    public void testPhaseCombining()
    {
        /*
        12:34:46.465 [main] [INFO ] 0: range(140490106 - 140500154) len(10048) phases(1000) preGene(1000)
        12:34:46.465 [main] [INFO ] 1: range(140500106 - 140508066) len(7960) phases(1000) preGene(0)
        12:34:46.465 [main] [INFO ] 2: range(140500137 - 140509812) len(9675) phases(1) preGene(0)
        12:34:46.465 [main] [INFO ] 3: range(140508067 - 140509010) len(943) phases(100) preGene(0)
        12:34:46.465 [main] [INFO ] 4: range(140509011 - 140509525) len(514) phases(1000) preGene(0)
        */

        List<GenePhaseRegion> simplePhases = Lists.newArrayList();
        simplePhases.add(new GenePhaseRegion("1", 140490106, 140500154, PHASE_1));
        simplePhases.get(0).setPreGene(true, PHASE_1);
        simplePhases.add(new GenePhaseRegion("1", 140490106, 140508066, PHASE_1));
        simplePhases.add(new GenePhaseRegion("1", 140500137, 140509812, PHASE_NON_CODING));
        simplePhases.add(new GenePhaseRegion("1", 140508067, 140509010, PHASE_0));
        simplePhases.add(new GenePhaseRegion("1", 140509011, 140509525, PHASE_1));

        List<GenePhaseRegion> combinedPhases = Lists.newArrayList();

        simplePhases.stream().forEach(x -> checkAddCombinedGenePhaseRegion(x, combinedPhases));

        assertFalse(validateSimpleVsCombinedPhaseRegions("1", simplePhases, combinedPhases));


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
    
    @Test
    public void testIntegratedCounts()
    {
        // test multiple chromosomes together for the various types of counts

        SvGeneTranscriptCollection geneTransCache = new SvGeneTranscriptCollection();

        int geneIndex = 1;
        String geneId, geneName;

        int transId = 1;
        long geneStart, geneEnd;

        List<EnsemblGeneData> geneList = Lists.newArrayList();
        List<TranscriptExonData> transExonList = null;

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

                transExonList = Lists.newArrayList();
                geneId = geneName = String.format("ESNG%04d", geneIndex++);
                createTransExons(transExonList, geneId, transId++, strand, geneExons1, new int[] {-1, 0, 1, 0, -1}, 100);
                addTransExonData(geneTransCache, geneId, transExonList);
                geneStart = transExonList.get(0).TransStart;
                geneEnd = transExonList.get(transExonList.size() - 1).TransEnd;
                geneList.add(createEnsemblGeneData(geneId, geneName, chromosome, strand, geneStart, geneEnd));

                transExonList = Lists.newArrayList();
                geneId = geneName = String.format("ESNG%04d", geneIndex++);
                createTransExons(transExonList, geneId, transId++, strand, geneExons2, new int[] {-1, 0, -1}, 100);
                addTransExonData(geneTransCache, geneId, transExonList);
                geneStart = transExonList.get(0).TransStart;
                geneEnd = transExonList.get(transExonList.size() - 1).TransEnd;
                geneList.add(createEnsemblGeneData(geneId, geneName, chromosome, strand, geneStart, geneEnd));

                // another transcript outside the DEL and DUP bucket length ranges
                transExonList = Lists.newArrayList();
                geneId = geneName = String.format("ESNG%04d", geneIndex++);
                createTransExons(transExonList, geneId, transId++, strand, geneExons3, new int[] {-1, 0, -1}, 100);
                addTransExonData(geneTransCache, geneId, transExonList);
                geneStart = transExonList.get(0).TransStart;
                geneEnd = transExonList.get(transExonList.size() - 1).TransEnd;
                geneList.add(createEnsemblGeneData(geneId, geneName, chromosome, strand, geneStart, geneEnd));

            }

            addGeneData(geneTransCache, chromosome, geneList);
            geneList = Lists.newArrayList();
        }

        // DELs can link sae gene and both closer genes, DUPs on the closer genes
        List<Long> delDupLengths = Lists.newArrayList((long)50, (long)5000, (long)50000);
        int shortInv = 20000;


        CohortExpFusions likelihoodCalc = new CohortExpFusions();
        likelihoodCalc.initialise(delDupLengths, shortInv, false);
        likelihoodCalc.generateGenePhasingCounts(geneTransCache, Lists.newArrayList(), Lists.newArrayList());
        likelihoodCalc.generateProximateFusionCounts();
        likelihoodCalc.generateNonProximateCounts();

        // expected overlap counts
        int intronLength = 1000;
        long shortDelOverlap = intronLength * intronLength;
        long medDelOverlap = 2 * intronLength * intronLength;

        final Map<String, List<GeneRangeData>> chrGeneDataMap = likelihoodCalc.getChrGeneRangeDataMap();
        assertEquals(3, chrGeneDataMap.size());

        List<GeneRangeData> geneRangeList = chrGeneDataMap.get("1");
        assertEquals(6, geneRangeList.size());

        GeneRangeData geneData = geneRangeList.get(0);
        assertEquals(4, geneData.getCombinedPhaseRegions().size());
        assertEquals(4, geneData.getPhaseRegions().size());

        assertEquals(2, geneData.getDelFusionBaseCounts().size());

        assertTrue(approxEqual(shortDelOverlap, geneData.getDelFusionBaseCounts().get(0), 0.01));

    }

    private boolean approxEqual(long bases1, long bases2, double allowPerc)
    {
        if((bases1 == 0) != (bases2 == 0))
            return false;

        double diff1 = (abs(bases1 - bases2) / (double)bases1);
        double diff2 = (abs(bases1 - bases2) / (double)bases2);

        return diff1 <= allowPerc && diff2 <= allowPerc;
    }

}
