package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.isofox.TestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.isofox.TestUtils.POS_STRAND;
import static com.hartwig.hmftools.isofox.TestUtils.createIsofoxConfig;
import static com.hartwig.hmftools.isofox.common.FragmentMatchType.LONG;
import static com.hartwig.hmftools.isofox.common.FragmentMatchType.SHORT;
import static com.hartwig.hmftools.isofox.common.FragmentMatchType.SPLICED;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.common.FragmentMatchType.UNSPLICED;
import static com.hartwig.hmftools.isofox.expression.ExpectedRatesGenerator.createTransComboDataMap;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.isofox.adjusts.FragmentSize;
import com.hartwig.hmftools.isofox.common.FragmentMatchType;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.expression.CategoryCountsData;
import com.hartwig.hmftools.common.sigs.ExpectationMaxFit;
import com.hartwig.hmftools.isofox.expression.ExpectedRatesData;
import com.hartwig.hmftools.isofox.expression.ExpectedRatesGenerator;
import com.hartwig.hmftools.common.utils.Matrix;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.junit.Test;

public class ExpectedRatesTest
{
    @Test
    public void testRegionMatching()
    {
        IsofoxConfig config = createIsofoxConfig();
        int fragmentLength = 100;
        config.ReadLength = 20;

        ExpectedRatesGenerator expRatesCalc = ExpectedRatesGenerator.from(config);

        expRatesCalc.setFragmentLengthData(fragmentLength, 1);

        int transId1 = 1;
        TranscriptData transData = new TranscriptData(transId1, "TRANS01", GENE_NAME_1, true, (byte)1,
                0, 1000, null,null, "");

        transData.exons().add(new ExonData(transId1, 100, 200, 1, -1, -1));
        transData.exons().add(new ExonData(transId1, 300, 400, 2, -1, -1));
        transData.exons().add(new ExonData(transId1, 440, 449, 3, -1, -1));
        transData.exons().add(new ExonData(transId1, 460, 469, 4, -1, -1));
        transData.exons().add(new ExonData(transId1, 600, 800, 5, -1, -1));

        // fully contained fragment
        int startPos = 100;
        List<int[]> readRegions = Lists.newArrayList();
        List<int[]> spliceJunctions = Lists.newArrayList();
        FragmentMatchType matchType = expRatesCalc.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertEquals(SHORT, matchType);
        assertEquals(2, readRegions.size());
        assertEquals(100, readRegions.get(0)[SE_START]);
        assertEquals(119, readRegions.get(0)[SE_END]);
        assertEquals(180, readRegions.get(1)[SE_START]);
        assertEquals(199, readRegions.get(1)[SE_END]);
        assertTrue(expRatesCalc.readsSupportTranscript(transData, readRegions, matchType, spliceJunctions));

        // test for another transcript
        int transId2 = 2;
        TranscriptData transData2 = new TranscriptData(transId2, "TRANS02", GENE_NAME_1, true, POS_STRAND,
                0, 1000, null,null, "");

        transData2.exons().add(new ExonData(transId2, 90, 210, 1, -1, -1));
        transData2.exons().add(new ExonData(transId2, 300, 400, 2, -1, -1));
        assertTrue(expRatesCalc.readsSupportTranscript(transData2, readRegions, matchType, spliceJunctions));

        transData2.exons().clear();
        transData2.exons().add(new ExonData(transId2, 150, 250, 1, -1, -1));
        transData2.exons().add(new ExonData(transId2, 300, 400, 2, -1, -1));
        assertFalse(expRatesCalc.readsSupportTranscript(transData2, readRegions, matchType, spliceJunctions));

        // 2 fully contained reads but in 2 exons
        startPos = 371;
        matchType = expRatesCalc.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertEquals(LONG, matchType);
        assertEquals(2, readRegions.size());
        assertEquals(371, readRegions.get(0)[SE_START]);
        assertEquals(390, readRegions.get(0)[SE_END]);
        assertEquals(631, readRegions.get(1)[SE_START]);
        assertEquals(650, readRegions.get(1)[SE_END]);
        assertTrue(expRatesCalc.readsSupportTranscript(transData, readRegions, matchType, spliceJunctions));

        transData2.exons().clear();
        transData2.exons().add(new ExonData(transId2, 250, 400, 1, -1, -1));
        transData2.exons().add(new ExonData(transId2, 450, 550, 2, -1, -1));
        transData2.exons().add(new ExonData(transId2, 600, 700, 3, -1, -1));
        assertTrue(expRatesCalc.readsSupportTranscript(transData2, readRegions, matchType, spliceJunctions));

        transData2.exons().clear();
        transData2.exons().add(new ExonData(transId2, 250, 300, 1, -1, -1));
        transData2.exons().add(new ExonData(transId2, 350, 700, 2, -1, -1));
        assertTrue(expRatesCalc.readsSupportTranscript(transData2, readRegions, matchType, spliceJunctions));

        // within an exon then spanning a junction
        startPos = 150;
        fragmentLength = 32 + 85 + 2 * config.ReadLength;
        expRatesCalc.setFragmentLengthData(fragmentLength, 1);

        matchType = expRatesCalc.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertEquals(SPLICED, matchType);
        assertEquals(3, readRegions.size());
        assertEquals(150, readRegions.get(0)[SE_START]);
        assertEquals(169, readRegions.get(0)[SE_END]);
        assertEquals(387, readRegions.get(1)[SE_START]);
        assertEquals(400, readRegions.get(1)[SE_END]);
        assertEquals(440, readRegions.get(2)[SE_START]);
        assertEquals(445, readRegions.get(2)[SE_END]);
        assertTrue(expRatesCalc.readsSupportTranscript(transData, readRegions, matchType, spliceJunctions));

        // start spanning a junction, ends within an exon
        startPos = 191;
        fragmentLength = 40 + 2 * config.ReadLength;
        expRatesCalc.setFragmentLengthData(fragmentLength, 1);

        matchType = expRatesCalc.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertEquals(SPLICED, matchType);
        assertEquals(3, readRegions.size());
        assertEquals(191, readRegions.get(0)[SE_START]);
        assertEquals(200, readRegions.get(0)[SE_END]);
        assertEquals(300, readRegions.get(1)[SE_START]);
        assertEquals(309, readRegions.get(1)[SE_END]);
        assertEquals(350, readRegions.get(2)[SE_START]);
        assertEquals(369, readRegions.get(2)[SE_END]);
        assertTrue(expRatesCalc.readsSupportTranscript(transData, readRegions, matchType, spliceJunctions));

        transData2.exons().clear();
        transData2.exons().add(new ExonData(transId2, 150, 200, 1, -1, -1));
        transData2.exons().add(new ExonData(transId2, 300, 320, 2, -1, -1));
        transData2.exons().add(new ExonData(transId2, 340, 380, 3, -1, -1));
        assertTrue(expRatesCalc.readsSupportTranscript(transData2, readRegions, matchType, spliceJunctions));

        // invalid since an exon is skipped in the splicing read
        transData2.exons().clear();
        transData2.exons().add(new ExonData(transId2, 150, 200, 1, -1, -1));
        transData2.exons().add(new ExonData(transId2, 210, 290, 2, -1, -1));
        transData2.exons().add(new ExonData(transId2, 300, 450, 3, -1, -1));
        assertFalse(expRatesCalc.readsSupportTranscript(transData2, readRegions, matchType, spliceJunctions));

        // 2 sets of exon junctions
        startPos = 191;
        fragmentLength = 81 + 5 + 2 * config.ReadLength;
        expRatesCalc.setFragmentLengthData(fragmentLength, 1);

        matchType = expRatesCalc.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertEquals(SPLICED, matchType);
        assertEquals(5, readRegions.size());
        assertEquals(191, readRegions.get(0)[SE_START]);
        assertEquals(200, readRegions.get(0)[SE_END]);
        assertEquals(300, readRegions.get(1)[SE_START]);
        assertEquals(309, readRegions.get(1)[SE_END]);
        assertEquals(396, readRegions.get(2)[SE_START]);
        assertEquals(400, readRegions.get(2)[SE_END]);
        assertEquals(440, readRegions.get(3)[SE_START]);
        assertEquals(449, readRegions.get(3)[SE_END]);
        assertEquals(460, readRegions.get(4)[SE_START]);
        assertEquals(464, readRegions.get(4)[SE_END]);
        assertTrue(expRatesCalc.readsSupportTranscript(transData, readRegions, matchType, spliceJunctions));

        // test cannot generate fragments past the end of the transcript
        startPos = 465;
        expRatesCalc.setFragmentLengthData(250, 1);
        matchType = expRatesCalc.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertTrue(readRegions.isEmpty());

        // test fragment sizes less than 2 or even 1 read length
        fragmentLength = 30;
        expRatesCalc.setFragmentLengthData(fragmentLength, 1);

        startPos = 150;
        matchType = expRatesCalc.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertEquals(SHORT, matchType);
        assertEquals(1, readRegions.size());
        assertEquals(150, readRegions.get(0)[SE_START]);
        assertEquals(179, readRegions.get(0)[SE_END]);
        assertTrue(expRatesCalc.readsSupportTranscript(transData, readRegions, matchType, spliceJunctions));
        assertTrue(expRatesCalc.readsSupportTranscript(transData2, readRegions, matchType, spliceJunctions));

        startPos = 190;
        matchType = expRatesCalc.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertEquals(SPLICED, matchType);
        assertEquals(2, readRegions.size());
        assertEquals(190, readRegions.get(0)[SE_START]);
        assertEquals(200, readRegions.get(0)[SE_END]);
        assertEquals(300, readRegions.get(1)[SE_START]);
        assertEquals(318, readRegions.get(1)[SE_END]);
        assertTrue(expRatesCalc.readsSupportTranscript(transData, readRegions, matchType, spliceJunctions));

        transData2.exons().clear();
        transData2.exons().add(new ExonData(transId2, 150, 200, 1, -1, -1));
        transData2.exons().add(new ExonData(transId2, 300, 320, 2, -1, -1));
        transData2.exons().add(new ExonData(transId2, 340, 380, 3, -1, -1));
        assertTrue(expRatesCalc.readsSupportTranscript(transData2, readRegions, matchType, spliceJunctions));

        // now with a fragment length less than the read length
        fragmentLength = 15;
        expRatesCalc.setFragmentLengthData(fragmentLength, 1);

        startPos = 185;
        matchType = expRatesCalc.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertEquals(SHORT, matchType);
        assertEquals(1, readRegions.size());
        assertEquals(185, readRegions.get(0)[SE_START]);
        assertEquals(199, readRegions.get(0)[SE_END]);
        assertTrue(expRatesCalc.readsSupportTranscript(transData, readRegions, matchType, spliceJunctions));
        assertTrue(expRatesCalc.readsSupportTranscript(transData2, readRegions, matchType, spliceJunctions));

        startPos = 195;
        matchType = expRatesCalc.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertEquals(SPLICED, matchType);
        assertEquals(2, readRegions.size());
        assertEquals(195, readRegions.get(0)[SE_START]);
        assertEquals(200, readRegions.get(0)[SE_END]);
        assertEquals(300, readRegions.get(1)[SE_START]);
        assertEquals(308, readRegions.get(1)[SE_END]);
        assertTrue(expRatesCalc.readsSupportTranscript(transData, readRegions, matchType, spliceJunctions));

        // test with a transcript shorter than the read length
        fragmentLength = 50;
        config.ReadLength = 150;
        expRatesCalc.setFragmentLengthData(fragmentLength, 1);

        int transId3 = 3;
        TranscriptData transData3 = new TranscriptData(transId3, "TRANS03", GENE_NAME_1, true, POS_STRAND,
                0, 1000, null,null, "");

        transData3.exons().add(new ExonData(transId3, 100, 239, 1, -1, -1));

        startPos = 100;
        matchType = expRatesCalc.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertEquals(SHORT, matchType);
        assertEquals(1, readRegions.size());
        assertEquals(100, readRegions.get(0)[SE_START]);
        assertEquals(149, readRegions.get(0)[SE_END]);
        assertTrue(expRatesCalc.readsSupportTranscript(transData3, readRegions, matchType, spliceJunctions));

        // and with a fragment length longer than the transcript
        fragmentLength = 150;
        config.ReadLength = 150;
        expRatesCalc.setFragmentLengthData(fragmentLength, 1);

        startPos = 100;
        matchType = expRatesCalc.generateImpliedFragment(transData3, startPos, readRegions, spliceJunctions);

        assertEquals(UNSPLICED, matchType);
        assertEquals(0, readRegions.size());
        assertFalse(expRatesCalc.readsSupportTranscript(transData3, readRegions, matchType, spliceJunctions));
    }

    @Test
    public void testSingleTranscriptCounts()
    {
        IsofoxConfig config = createIsofoxConfig();
        config.ReadLength = 10;
        config.FragmentSizeData.add(new FragmentSize(30, 1));

        ExpectedRatesGenerator expRatesCalc = ExpectedRatesGenerator.from(config);

        String geneId = "GENE01";

        GeneData geneData = new GeneData(geneId, geneId, "1", POS_STRAND, 100, 1000, "");

        int transId = 1;
        String transName = "TRANS01";

        TranscriptData transData = new TranscriptData(transId, transName, geneId, true, POS_STRAND,
                100, 414, null,null, "");

        transData.exons().add(new ExonData(transId, 100, 158, 1, -1, -1));
        transData.exons().add(new ExonData(transId, 228, 286, 2, -1, -1));
        transData.exons().add(new ExonData(transId, 356, 414, 3, -1, -1));

        GeneReadData geneReadData = new GeneReadData(geneData);

        List<TranscriptData> transcripts = Lists.newArrayList(transData);

        geneReadData.setTranscripts(transcripts);

        GeneCollection genes = new GeneCollection(0, Lists.newArrayList(geneReadData));

        expRatesCalc.generateExpectedRates(genes);

        Map<String,List<CategoryCountsData>> transComboData = createTransComboDataMap(expRatesCalc.getTransComboData());
        assertEquals(2, transComboData.size());
        String transIdStr = String.valueOf(transId);
        assertTrue(transComboData.containsKey(transIdStr));
        assertTrue(transComboData.containsKey(geneId));

        List<CategoryCountsData> tcDataList = transComboData.get(transIdStr);
        assertEquals(2, tcDataList.size());

        List<Integer> tranIds = Lists.newArrayList(transId);
        List<String> unsplicedGenes = Lists.newArrayList(geneId);

        CategoryCountsData tcData = findMatchingData(tranIds, unsplicedGenes, tcDataList);
        assertTrue(tcData != null);
        assertEquals(90, tcData.fragmentCount(), 0.01);

        tcData = findMatchingData(tranIds, Lists.newArrayList(), tcDataList);
        assertEquals(58, tcData.fragmentCount(), 0.01);

        tcDataList = transComboData.get(geneId);
        assertEquals(2, tcDataList.size());

        tcData = findMatchingData(tranIds, unsplicedGenes, tcDataList);
        assertTrue(tcData != null);
        assertEquals(90, tcData.fragmentCount(), 0.01);

        tranIds.clear();
        tcData = findMatchingData(tranIds, unsplicedGenes, tcDataList);
        assertTrue(tcData != null);
        assertEquals(196, tcData.fragmentCount(), 0.01);
    }

    private static CategoryCountsData findMatchingData(
            final List<Integer> transcripts, final List<String> unsplicedGenes, final List<CategoryCountsData> dataList)
    {
        return dataList.stream().filter(x -> x.matches(transcripts, unsplicedGenes)).findFirst().orElse(null);
    }

    @Test
    public void testMultipleTranscriptCounts()
    {
        IsofoxConfig config = createIsofoxConfig();
        config.FragmentSizeData.add(new FragmentSize(30, 1));
        config.ReadLength = 10;

        ExpectedRatesGenerator expRatesCalc = ExpectedRatesGenerator.from(config);

        String geneId = "GENE01";

        GeneData geneData = new GeneData(geneId, geneId, "1", POS_STRAND, 100, 1000, "");

        int transId1 = 1;
        String transName1 = "TRANS01";

        TranscriptData transData1 = new TranscriptData(transId1, transName1, geneId, true, POS_STRAND,
                100, 600, null,null, "");

        transData1.exons().add(new ExonData(transId1, 100, 200, 1, -1, -1));
        transData1.exons().add(new ExonData(transId1, 300, 400, 2, -1, -1));
        transData1.exons().add(new ExonData(transId1, 500, 600, 3, -1, -1));

        int transId2 = 2;
        String transName2 = "TRANS02";

        TranscriptData transData2 = new TranscriptData(transId2, transName2, geneId, true, POS_STRAND,
                150, 1000, null,null, "");

        transData2.exons().add(new ExonData(transId2, 150, 200, 1, -1, -1));
        transData2.exons().add(new ExonData(transId2, 300, 450, 2, -1, -1));
        transData2.exons().add(new ExonData(transId2, 700, 800, 3, -1, -1));
        transData2.exons().add(new ExonData(transId2, 900, 1000, 4, -1, -1));

        GeneReadData geneReadData = new GeneReadData(geneData);

        List<TranscriptData> transcripts = Lists.newArrayList(transData1, transData2);

        geneReadData.setTranscripts(transcripts);

        GeneCollection genes = new GeneCollection(0, Lists.newArrayList(geneReadData));

        List<int[]> commonExonicRegions = genes.getCommonExonicRegions();
        assertEquals(5, commonExonicRegions.size());
        assertEquals(100, commonExonicRegions.get(0)[SE_START]);
        assertEquals(200, commonExonicRegions.get(0)[SE_END]);
        assertEquals(300, commonExonicRegions.get(1)[SE_START]);
        assertEquals(450, commonExonicRegions.get(1)[SE_END]);
        assertEquals(500, commonExonicRegions.get(2)[SE_START]);
        assertEquals(600, commonExonicRegions.get(2)[SE_END]);
        assertEquals(700, commonExonicRegions.get(3)[SE_START]);
        assertEquals(800, commonExonicRegions.get(3)[SE_END]);
        assertEquals(900, commonExonicRegions.get(4)[SE_START]);
        assertEquals(1000, commonExonicRegions.get(4)[SE_END]);

        expRatesCalc.generateExpectedRates(genes);

        Map<String,List<CategoryCountsData>> transComboData = expRatesCalc.getTransComboDataMap();
        assertEquals(3, transComboData.size());
        assertTrue(transComboData.containsKey(transName1));
        assertTrue(transComboData.containsKey(transName2));
        assertTrue(transComboData.containsKey(geneId));

        List<CategoryCountsData> tcDataList = transComboData.get(transName1);
        assertEquals(4, tcDataList.size());

        List<Integer> tranIds = Lists.newArrayList(transId1);
        List<String> unsplicedGenes = Lists.newArrayList(geneId);

        // short trans 1 or unspliced
        CategoryCountsData tcData = findMatchingData(tranIds, unsplicedGenes, tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.fragmentCount() > 0);

        // only spliced trans 1
        tcData = findMatchingData(tranIds, Lists.newArrayList(), tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.fragmentCount() > 0);

        // short trans 1 or 2 or unspliced
        tranIds = Lists.newArrayList(transId1, transId2);
        tcData = findMatchingData(tranIds, unsplicedGenes, tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.fragmentCount() > 0);

        // trans 1 or 2 spliced
        tranIds = Lists.newArrayList(transId1, transId2);
        tcData = findMatchingData(tranIds, Lists.newArrayList(), tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.fragmentCount() > 0);

        tcDataList = transComboData.get(transName2);
        assertEquals(4, tcDataList.size());

        tranIds = Lists.newArrayList(transId2);

        tcData = findMatchingData(tranIds, unsplicedGenes, tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.fragmentCount() > 0);

        tranIds = Lists.newArrayList(transId1, transId2);
        tcData = findMatchingData(tranIds, unsplicedGenes, tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.fragmentCount() > 0);

        tcDataList = transComboData.get(geneId);
        assertEquals(4, tcDataList.size());

        tranIds.clear();
        tcData = findMatchingData(tranIds, unsplicedGenes, tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.fragmentCount() > 0);
    }

    @Test
    public void testSingleExonicRegions()
    {
        IsofoxConfig config = createIsofoxConfig();
        config.FragmentSizeData.add(new FragmentSize(30, 1));
        config.ReadLength = 10;

        ExpectedRatesGenerator expRatesCalc = ExpectedRatesGenerator.from(config);

        String geneId = "GENE01";

        GeneData geneData = new GeneData(geneId, geneId, "1", POS_STRAND, 100, 400, "");

        // first the single transcript and single exon
        int transId = 1;
        String transName1 = "TRANS01";

        TranscriptData transData1 = new TranscriptData(transId, transName1, geneId, true, POS_STRAND,
                100, 300, null,null, "");

        transData1.exons().add(new ExonData(transId, 100, 300, 1, -1, -1));

        GeneReadData geneReadData = new GeneReadData(geneData);

        List<TranscriptData> transcripts = Lists.newArrayList(transData1);

        geneReadData.setTranscripts(transcripts);

        GeneCollection genes = new GeneCollection(0, Lists.newArrayList(geneReadData));

        expRatesCalc.generateExpectedRates(genes);

        Map<String,List<CategoryCountsData>> transComboData = expRatesCalc.getTransComboDataMap();
        assertEquals(2, transComboData.size());
        assertTrue(transComboData.containsKey(transName1));
        assertTrue(transComboData.containsKey(geneId));

        ExpectedRatesData erData = expRatesCalc.getExpectedRatesData();

        Matrix rates = erData.getTranscriptDefinitions();
        assertEquals(2, rates.Cols);
        assertEquals(2, rates.Rows);

        assertEquals(2, erData.Categories.size());
        assertEquals(geneId, erData.Categories.get(1));

        final String transCategory = transComboData.get(transName1).get(0).combinedKey();

        int transCatIndex = erData.getCategoryIndex(transCategory);
        final String transIdStr = String.valueOf(transId);
        int transIndex = erData.getTranscriptIndex(transIdStr);
        assertEquals(transCatIndex, 0);
        assertEquals(transIndex, 0);

        assertEquals(1, rates.get(transCatIndex, transIndex), 0.001); // 100% in the actual trans' count
        assertEquals(0, rates.get(1, 0), 0.001);
        assertEquals(0, rates.get(0, 1), 0.001);
        assertEquals(0, rates.get(1, 1), 0.001);

        int transId2 = 2;
        String transName2 = "TRANS02";

        TranscriptData transData2 = new TranscriptData(transId2, transName2, geneId, true, POS_STRAND,
                200, 400, null,null, "");

        transData2.exons().add(new ExonData(transId, 200, 400, 1, -1, -1));

        geneReadData = new GeneReadData(geneData);

        transcripts = Lists.newArrayList(transData1, transData2);

        geneReadData.setTranscripts(transcripts);

        genes = new GeneCollection(0, Lists.newArrayList(geneReadData));

        List<int[]> commonExonicRegions = genes.getCommonExonicRegions();
        assertEquals(1, commonExonicRegions.size());
        assertEquals(100, commonExonicRegions.get(0)[SE_START]);
        assertEquals(400, commonExonicRegions.get(0)[SE_END]);

        expRatesCalc.generateExpectedRates(genes);

        transComboData = expRatesCalc.getTransComboDataMap();
        assertEquals(3, transComboData.size());
        assertTrue(transComboData.containsKey(transName1));
        assertTrue(transComboData.containsKey(transName2));
        assertTrue(transComboData.containsKey(geneId));
    }

    @Test
    public void testOverlappingGeneCounts()
    {
        IsofoxConfig config = createIsofoxConfig();
        config.FragmentSizeData.add(new FragmentSize(30, 1));
        config.ReadLength = 10;

        ExpectedRatesGenerator expRatesCalc = ExpectedRatesGenerator.from(config);

        String geneId1 = "GENE01";
        GeneData geneData1 = new GeneData(geneId1, geneId1, "1", (byte) 1, 100, 600, "");

        int transId1 = 1;
        String transName1 = "TRANS01";

        TranscriptData transData1 = new TranscriptData(transId1, transName1, geneId1, true, (byte) 1,
                100, 600, null, null, "");

        transData1.exons().add(new ExonData(transId1, 100, 200, 1, -1, -1));
        transData1.exons().add(new ExonData(transId1, 300, 400, 2, -1, -1));
        transData1.exons().add(new ExonData(transId1, 500, 600, 3, -1, -1));

        GeneReadData geneReadData1 = new GeneReadData(geneData1);
        geneReadData1.setTranscripts(Lists.newArrayList(transData1));

        String geneId2 = "GENE02";
        GeneData geneData2 = new GeneData(geneId2, geneId2, "1", (byte) 1, 150, 1000, "");

        int transId2 = 2;
        String transName2 = "TRANS02";

        TranscriptData transData2 = new TranscriptData(transId2, transName2, geneId2, true, (byte) 1,
                150, 1000, null, null, "");

        transData2.exons().add(new ExonData(transId2, 150, 200, 1, -1, -1));
        transData2.exons().add(new ExonData(transId2, 300, 450, 2, -1, -1));
        transData2.exons().add(new ExonData(transId2, 700, 800, 3, -1, -1));
        transData2.exons().add(new ExonData(transId2, 900, 1000, 4, -1, -1));

        GeneReadData geneReadData2 = new GeneReadData(geneData2);

        geneReadData2.setTranscripts(Lists.newArrayList(transData2));

        GeneCollection genes = new GeneCollection(0, Lists.newArrayList(geneReadData1, geneReadData2));

        List<int[]> commonExonicRegions = genes.getCommonExonicRegions();
        assertEquals(5, commonExonicRegions.size());

        expRatesCalc.generateExpectedRates(genes);

        Map<String,List<CategoryCountsData>> transComboData = expRatesCalc.getTransComboDataMap();
        assertEquals(4, transComboData.size());
        assertTrue(transComboData.containsKey(transName1));
        assertTrue(transComboData.containsKey(transName2));
        assertTrue(transComboData.containsKey(geneId1));
        assertTrue(transComboData.containsKey(geneId2));

        List<CategoryCountsData> tcDataList = transComboData.get(transName1);
        assertEquals(5, tcDataList.size());

        // trans 1 by itself
        CategoryCountsData tcData = findMatchingData(Lists.newArrayList(transId1), Lists.newArrayList(), tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.fragmentCount() > 0);

        tcData = findMatchingData(Lists.newArrayList(transId1, transId2), Lists.newArrayList(), tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.fragmentCount() > 0);

        tcData = findMatchingData(Lists.newArrayList(transId1), Lists.newArrayList(geneId1), tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.fragmentCount() > 0);

        tcData = findMatchingData(Lists.newArrayList(transId1, transId2), Lists.newArrayList(geneId1, geneId2), tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.fragmentCount() > 0);

        tcData = findMatchingData(Lists.newArrayList(transId1), Lists.newArrayList(geneId1, geneId2), tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.fragmentCount() > 0);

        tcDataList = transComboData.get(transName2);
        assertEquals(5, tcDataList.size());

        tcDataList = transComboData.get(geneId1);
        assertEquals(5, tcDataList.size());

        tcData = findMatchingData(Lists.newArrayList(), Lists.newArrayList(geneId1), tcDataList);
        assertTrue(tcData == null);

        tcData = findMatchingData(Lists.newArrayList(transId1), Lists.newArrayList(geneId1), tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.fragmentCount() > 0);

        tcData = findMatchingData(Lists.newArrayList(transId1), Lists.newArrayList(geneId1, geneId2), tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.fragmentCount() > 0);

        tcData = findMatchingData(Lists.newArrayList(), Lists.newArrayList(geneId1, geneId2), tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.fragmentCount() > 0);

        tcData = findMatchingData(Lists.newArrayList(transId1, transId2), Lists.newArrayList(geneId1, geneId2), tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.fragmentCount() > 0);

        tcData = findMatchingData(Lists.newArrayList(transId2), Lists.newArrayList(geneId1, geneId2), tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.fragmentCount() > 0);

        tcDataList = transComboData.get(geneId2);
        assertEquals(6, tcDataList.size());

        tcData = findMatchingData(Lists.newArrayList(), Lists.newArrayList(geneId2), tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.fragmentCount() > 0);
    }

    @Test
    public void testExpectationMaxFit()
    {
        Configurator.setRootLevel(Level.DEBUG);

        int categoryCount = 3;
        int transCount = 2;

        Matrix sigs = new Matrix(categoryCount, transCount);
        double[] transSig1 = {0.2, 0.8, 0};
        double[] transSig2 = {0.4, 0, 0.6};

        sigs.setCol(0, transSig1);
        sigs.setCol(1, transSig2);

        double[] transCounts = new double[categoryCount];

        transCounts[0] = 5;
        transCounts[1] = 4;
        transCounts[2] = 6;

        double[] allocations = ExpectationMaxFit.performFit(transCounts, sigs);

        assertEquals(5.002, allocations[0], 0.001);
        assertEquals(9.998, allocations[1], 0.001);

        transCounts[0] = 5;
        transCounts[1] = 4;
        transCounts[2] = 7;

        allocations = ExpectationMaxFit.performFit(transCounts, sigs);

        assertEquals(4.905, allocations[0], 0.001);
        assertEquals(11.095, allocations[1], 0.001);

    }

}
