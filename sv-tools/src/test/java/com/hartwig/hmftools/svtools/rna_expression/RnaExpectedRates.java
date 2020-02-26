package com.hartwig.hmftools.svtools.rna_expression;

import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.ExpectedExpressionRates.UNSPLICED_ID;
import static com.hartwig.hmftools.svtools.rna_expression.FragmentMatchType.LONG;
import static com.hartwig.hmftools.svtools.rna_expression.FragmentMatchType.SHORT;
import static com.hartwig.hmftools.svtools.rna_expression.FragmentMatchType.SPLICED;
import static com.hartwig.hmftools.svtools.rna_expression.TranscriptComboData.findMatchingData;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.sig_analyser.common.SigMatrix;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configuration;
import org.apache.logging.log4j.core.config.Configurator;
import org.junit.Test;

public class RnaExpectedRates
{
    @Test
    public void testRegionMatching()
    {
        RnaExpConfig config = new RnaExpConfig();
        int fragmentLength = 100;
        config.ReadLength = 20;

        ExpectedExpressionRates eeRates = ExpectedExpressionRates.from(config);

        eeRates.setFragmentLengthData(fragmentLength, 1);

        int transId1 = 1;
        TranscriptData transData = new TranscriptData(transId1, "TRANS01", "GENE01", true, (byte)1,
                0, 1000, null,null, "");

        transData.exons().add(new ExonData(transId1, 100, 200, 1, -1, -1));
        transData.exons().add(new ExonData(transId1, 300, 400, 2, -1, -1));
        transData.exons().add(new ExonData(transId1, 440, 449, 3, -1, -1));
        transData.exons().add(new ExonData(transId1, 460, 469, 4, -1, -1));
        transData.exons().add(new ExonData(transId1, 600, 800, 5, -1, -1));

        // fully contained fragment
        long startPos = 100;
        List<long[]> readRegions = Lists.newArrayList();
        List<long[]> spliceJunctions = Lists.newArrayList();
        FragmentMatchType matchType = eeRates.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertEquals(SHORT, matchType);
        assertEquals(2, readRegions.size());
        assertEquals(100, readRegions.get(0)[SE_START]);
        assertEquals(119, readRegions.get(0)[SE_END]);
        assertEquals(180, readRegions.get(1)[SE_START]);
        assertEquals(199, readRegions.get(1)[SE_END]);
        assertTrue(eeRates.readsSupportFragment(transData, readRegions, matchType, spliceJunctions));

        // test for another transcript
        int transId2 = 2;
        TranscriptData transData2 = new TranscriptData(2, "TRANS02", "GENE01", true, (byte)1,
                0, 1000, null,null, "");

        transData2.exons().add(new ExonData(transId2, 90, 210, 1, -1, -1));
        transData2.exons().add(new ExonData(transId2, 300, 400, 2, -1, -1));
        assertTrue(eeRates.readsSupportFragment(transData2, readRegions, matchType, spliceJunctions));

        transData2.exons().clear();
        transData2.exons().add(new ExonData(transId2, 150, 250, 1, -1, -1));
        transData2.exons().add(new ExonData(transId2, 300, 400, 2, -1, -1));
        assertFalse(eeRates.readsSupportFragment(transData2, readRegions, matchType, spliceJunctions));

        // 2 fully contained reads but in 2 exons
        startPos = 371;
        matchType = eeRates.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertEquals(LONG, matchType);
        assertEquals(2, readRegions.size());
        assertEquals(371, readRegions.get(0)[SE_START]);
        assertEquals(390, readRegions.get(0)[SE_END]);
        assertEquals(631, readRegions.get(1)[SE_START]);
        assertEquals(650, readRegions.get(1)[SE_END]);
        assertTrue(eeRates.readsSupportFragment(transData, readRegions, matchType, spliceJunctions));

        transData2.exons().clear();
        transData2.exons().add(new ExonData(transId2, 250, 400, 1, -1, -1));
        transData2.exons().add(new ExonData(transId2, 450, 550, 2, -1, -1));
        transData2.exons().add(new ExonData(transId2, 600, 700, 3, -1, -1));
        assertTrue(eeRates.readsSupportFragment(transData2, readRegions, matchType, spliceJunctions));

        transData2.exons().clear();
        transData2.exons().add(new ExonData(transId2, 250, 300, 1, -1, -1));
        transData2.exons().add(new ExonData(transId2, 350, 700, 2, -1, -1));
        assertTrue(eeRates.readsSupportFragment(transData2, readRegions, matchType, spliceJunctions));

        // within an exon then spanning a junction
        startPos = 150;
        fragmentLength = 32 + 85 + 2 * config.ReadLength;
        eeRates.setFragmentLengthData(fragmentLength, 1);

        matchType = eeRates.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertEquals(SPLICED, matchType);
        assertEquals(3, readRegions.size());
        assertEquals(150, readRegions.get(0)[SE_START]);
        assertEquals(169, readRegions.get(0)[SE_END]);
        assertEquals(387, readRegions.get(1)[SE_START]);
        assertEquals(400, readRegions.get(1)[SE_END]);
        assertEquals(440, readRegions.get(2)[SE_START]);
        assertEquals(445, readRegions.get(2)[SE_END]);
        assertTrue(eeRates.readsSupportFragment(transData, readRegions, matchType, spliceJunctions));

        // start spanning a junction, ends within an exon
        startPos = 191;
        fragmentLength = 40 + 2 * config.ReadLength;
        eeRates.setFragmentLengthData(fragmentLength, 1);

        matchType = eeRates.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertEquals(SPLICED, matchType);
        assertEquals(3, readRegions.size());
        assertEquals(191, readRegions.get(0)[SE_START]);
        assertEquals(200, readRegions.get(0)[SE_END]);
        assertEquals(300, readRegions.get(1)[SE_START]);
        assertEquals(309, readRegions.get(1)[SE_END]);
        assertEquals(350, readRegions.get(2)[SE_START]);
        assertEquals(369, readRegions.get(2)[SE_END]);
        assertTrue(eeRates.readsSupportFragment(transData, readRegions, matchType, spliceJunctions));

        transData2.exons().clear();
        transData2.exons().add(new ExonData(transId2, 150, 200, 1, -1, -1));
        transData2.exons().add(new ExonData(transId2, 300, 320, 2, -1, -1));
        transData2.exons().add(new ExonData(transId2, 340, 380, 3, -1, -1));
        assertTrue(eeRates.readsSupportFragment(transData2, readRegions, matchType, spliceJunctions));

        // invalid since an exon is skipped in the splicing read
        transData2.exons().clear();
        transData2.exons().add(new ExonData(transId2, 150, 200, 1, -1, -1));
        transData2.exons().add(new ExonData(transId2, 210, 290, 2, -1, -1));
        transData2.exons().add(new ExonData(transId2, 300, 450, 3, -1, -1));
        assertFalse(eeRates.readsSupportFragment(transData2, readRegions, matchType, spliceJunctions));

        // 2 sets of exon junctions
        startPos = 191;
        fragmentLength = 81 + 5 + 2 * config.ReadLength;
        eeRates.setFragmentLengthData(fragmentLength, 1);

        matchType = eeRates.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

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
        assertTrue(eeRates.readsSupportFragment(transData, readRegions, matchType, spliceJunctions));

        // test cannot generate fragments past the end of the transcript
        startPos = 465;
        eeRates.setFragmentLengthData(250, 1);
        matchType = eeRates.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertTrue(readRegions.isEmpty());

        // test fragment sizes less than 2 or even 1 read length
        fragmentLength = 30;
        eeRates.setFragmentLengthData(fragmentLength, 1);

        startPos = 150;
        matchType = eeRates.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertEquals(SHORT, matchType);
        assertEquals(1, readRegions.size());
        assertEquals(150, readRegions.get(0)[SE_START]);
        assertEquals(179, readRegions.get(0)[SE_END]);
        assertTrue(eeRates.readsSupportFragment(transData, readRegions, matchType, spliceJunctions));
        assertTrue(eeRates.readsSupportFragment(transData2, readRegions, matchType, spliceJunctions));

        startPos = 190;
        matchType = eeRates.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertEquals(SPLICED, matchType);
        assertEquals(2, readRegions.size());
        assertEquals(190, readRegions.get(0)[SE_START]);
        assertEquals(200, readRegions.get(0)[SE_END]);
        assertEquals(300, readRegions.get(1)[SE_START]);
        assertEquals(318, readRegions.get(1)[SE_END]);
        assertTrue(eeRates.readsSupportFragment(transData, readRegions, matchType, spliceJunctions));

        transData2.exons().clear();
        transData2.exons().add(new ExonData(transId2, 150, 200, 1, -1, -1));
        transData2.exons().add(new ExonData(transId2, 300, 320, 2, -1, -1));
        transData2.exons().add(new ExonData(transId2, 340, 380, 3, -1, -1));
        assertTrue(eeRates.readsSupportFragment(transData2, readRegions, matchType, spliceJunctions));

        // now with a fragment length less than the read length
        fragmentLength = 15;
        eeRates.setFragmentLengthData(fragmentLength, 1);

        startPos = 185;
        matchType = eeRates.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertEquals(SHORT, matchType);
        assertEquals(1, readRegions.size());
        assertEquals(185, readRegions.get(0)[SE_START]);
        assertEquals(199, readRegions.get(0)[SE_END]);
        assertTrue(eeRates.readsSupportFragment(transData, readRegions, matchType, spliceJunctions));
        assertTrue(eeRates.readsSupportFragment(transData2, readRegions, matchType, spliceJunctions));

        startPos = 195;
        matchType = eeRates.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertEquals(SPLICED, matchType);
        assertEquals(2, readRegions.size());
        assertEquals(195, readRegions.get(0)[SE_START]);
        assertEquals(200, readRegions.get(0)[SE_END]);
        assertEquals(300, readRegions.get(1)[SE_START]);
        assertEquals(308, readRegions.get(1)[SE_END]);
        assertTrue(eeRates.readsSupportFragment(transData, readRegions, matchType, spliceJunctions));
    }

    @Test
    public void testSingleTranscriptCounts()
    {
        RnaExpConfig config = new RnaExpConfig();
        config.ReadLength = 10;
        config.ExpRateFragmentLengths.add(new int[] {30, 1});

        ExpectedExpressionRates eeRates = ExpectedExpressionRates.from(config);

        String geneId = "GENE01";

        EnsemblGeneData geneData = new EnsemblGeneData(geneId, geneId, "1", (byte)1, 100, 1000, "");

        int transId = 1;
        String transName = "TRANS01";

        TranscriptData transData = new TranscriptData(transId, transName, geneId, true, (byte)1,
                100, 414, null,null, "");

        transData.exons().add(new ExonData(transId, 100, 158, 1, -1, -1));
        transData.exons().add(new ExonData(transId, 228, 286, 2, -1, -1));
        transData.exons().add(new ExonData(transId, 356, 414, 3, -1, -1));

        GeneReadData geneReadData = new GeneReadData(geneData);

        List<TranscriptData> transcripts = Lists.newArrayList(transData);

        geneReadData.setTranscripts(transcripts);
        geneReadData.generateExonicRegions();

        eeRates.generateExpectedRates(geneReadData);

        Map<String,List<TranscriptComboData>> transComboData = eeRates.getTransComboData();
        assertEquals(2, transComboData.size());
        assertTrue(transComboData.containsKey(transId));
        assertTrue(transComboData.containsKey(UNSPLICED_ID));

        List<TranscriptComboData> tcDataList = transComboData.get(transId);
        assertEquals(1, tcDataList.size());

        List<Integer> tranIds = Lists.newArrayList(transId);

        TranscriptComboData tcData = findMatchingData(tranIds, tcDataList);
        assertTrue(tcData != null);
        assertEquals(90, tcData.getShortCount());
        assertEquals(58, tcData.getSplicedCount());

        tcDataList = transComboData.get(UNSPLICED_ID);
        assertEquals(2, tcDataList.size());

        tcData = findMatchingData(tranIds, tcDataList);
        assertTrue(tcData != null);
        assertEquals(87, tcData.getShortCount());
        assertEquals(0, tcData.getSplicedCount());

        tranIds.clear();
        tcData = findMatchingData(tranIds, tcDataList);
        assertTrue(tcData != null);
        assertEquals(196, tcData.getUnsplicedCount());
        assertEquals(tcData.getUnsplicedCount(), tcData.totalCount());
    }

    @Test
    public void testMultipleTranscriptCounts()
    {
        RnaExpConfig config = new RnaExpConfig();
        config.ExpRateFragmentLengths.add(new int[] {30, 1});
        config.ReadLength = 10;

        ExpectedExpressionRates eeRates = ExpectedExpressionRates.from(config);

        String geneId = "GENE01";

        EnsemblGeneData geneData = new EnsemblGeneData(geneId, geneId, "1", (byte)1, 100, 1000, "");

        int transId1 = 1;

        TranscriptData transData1 = new TranscriptData(1, "TRANS01", geneId, true, (byte)1,
                100, 600, null,null, "");

        transData1.exons().add(new ExonData(transId1, 100, 200, 1, -1, -1));
        transData1.exons().add(new ExonData(transId1, 300, 400, 2, -1, -1));
        transData1.exons().add(new ExonData(transId1, 500, 600, 3, -1, -1));

        int transId2 = 2;

        TranscriptData transData2 = new TranscriptData(1, "TRANS02", geneId, true, (byte)1,
                150, 1000, null,null, "");

        transData2.exons().add(new ExonData(transId1, 150, 200, 1, -1, -1));
        transData2.exons().add(new ExonData(transId1, 300, 450, 2, -1, -1));
        transData2.exons().add(new ExonData(transId1, 700, 800, 3, -1, -1));
        transData2.exons().add(new ExonData(transId1, 900, 1000, 4, -1, -1));

        GeneReadData geneReadData = new GeneReadData(geneData);

        List<TranscriptData> transcripts = Lists.newArrayList(transData1, transData2);

        geneReadData.setTranscripts(transcripts);
        geneReadData.generateExonicRegions();

        List<long[]> commonExonicRegions = geneReadData.getCommonExonicRegions();
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

        eeRates.generateExpectedRates(geneReadData);

        Map<String,List<TranscriptComboData>> transComboData = eeRates.getTransComboData();
        assertEquals(3, transComboData.size());
        assertTrue(transComboData.containsKey(transId1));
        assertTrue(transComboData.containsKey(transId2));
        assertTrue(transComboData.containsKey(UNSPLICED_ID));

        List<TranscriptComboData> tcDataList = transComboData.get(transId1);
        assertEquals(2, tcDataList.size());

        List<Integer> tranIds = Lists.newArrayList(transId1);

        TranscriptComboData tcData = findMatchingData(tranIds, tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.getShortCount() > 0);
        assertTrue(tcData.getSplicedCount() > 0);
        assertEquals(tcData.getSplicedCount() + tcData.getShortCount(), tcData.totalCount());

        tranIds = Lists.newArrayList(transId1, transId2);
        tcData = findMatchingData(tranIds, tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.getShortCount() > 0);
        assertTrue(tcData.getSplicedCount() > 0);

        tcDataList = transComboData.get(transId2);
        assertEquals(2, tcDataList.size());

        tranIds = Lists.newArrayList(transId2);

        tcData = findMatchingData(tranIds, tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.getShortCount() > 0);
        assertTrue(tcData.getSplicedCount() > 0);

        tranIds = Lists.newArrayList(transId1, transId2);
        tcData = findMatchingData(tranIds, tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.getShortCount() > 0);
        assertTrue(tcData.getSplicedCount() > 0);

        tcDataList = transComboData.get(UNSPLICED_ID);
        assertEquals(4, tcDataList.size());

        tranIds.clear();
        tcData = findMatchingData(tranIds, tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.getUnsplicedCount() > 0);
        assertEquals(tcData.getUnsplicedCount(), tcData.totalCount());
    }

    @Test
    public void testSingleExonicRegions()
    {
        RnaExpConfig config = new RnaExpConfig();
        config.ExpRateFragmentLengths.add(new int[] {30, 1});
        config.ReadLength = 10;

        ExpectedExpressionRates eeRates = ExpectedExpressionRates.from(config);

        String geneId = "GENE01";

        EnsemblGeneData geneData = new EnsemblGeneData(geneId, geneId, "1", (byte)1, 100, 400, "");

        // first the single transcript and single exon
        int transId = 1;

        TranscriptData transData1 = new TranscriptData(transId, "TRANS01", geneId, true, (byte)1,
                100, 300, null,null, "");

        transData1.exons().add(new ExonData(transId, 100, 300, 1, -1, -1));

        GeneReadData geneReadData = new GeneReadData(geneData);

        List<TranscriptData> transcripts = Lists.newArrayList(transData1);

        geneReadData.setTranscripts(transcripts);
        geneReadData.generateExonicRegions();
        eeRates.generateExpectedRates(geneReadData);

        Map<String,List<TranscriptComboData>> transComboData = eeRates.getTransComboData();
        assertEquals(2, transComboData.size());
        assertTrue(transComboData.containsKey(transId));
        assertTrue(transComboData.containsKey(UNSPLICED_ID));

        SigMatrix rates = eeRates.getTranscriptDefinitions();
        assertEquals(2, rates.Cols);
        assertEquals(2, rates.Rows);

        assertEquals(2, eeRates.getCategories().size());
        assertEquals(UNSPLICED_ID, eeRates.getCategories().get(0));

        assertEquals(0, rates.get(0, 0), 0.001);
        assertEquals(0, rates.get(1, 0), 0.001);
        assertEquals(0, rates.get(0, 1), 0.001);
        assertEquals(1, rates.get(1, 1), 0.001); // 100% in the actual trans' count

        String transId2 = "TRANS02";

        TranscriptData transData2 = new TranscriptData(1, transId2, geneId, true, (byte)1,
                200, 400, null,null, "");

        transData2.exons().add(new ExonData(transId, 200, 400, 1, -1, -1));

        geneReadData = new GeneReadData(geneData);

        transcripts = Lists.newArrayList(transData1, transData2);

        geneReadData.setTranscripts(transcripts);
        geneReadData.generateExonicRegions();

        List<long[]> commonExonicRegions = geneReadData.getCommonExonicRegions();
        assertEquals(1, commonExonicRegions.size());
        assertEquals(100, commonExonicRegions.get(0)[SE_START]);
        assertEquals(400, commonExonicRegions.get(0)[SE_END]);

        eeRates.generateExpectedRates(geneReadData);

        transComboData = eeRates.getTransComboData();
        assertEquals(3, transComboData.size());
        assertTrue(transComboData.containsKey(transId));
        assertTrue(transComboData.containsKey(transId2));
        assertTrue(transComboData.containsKey(UNSPLICED_ID));
    }

    @Test
    public void testExpectationMaxFit()
    {
        Configurator.setRootLevel(Level.DEBUG);

        int categoryCount = 3;
        int transCount = 2;

        SigMatrix sigs = new SigMatrix(categoryCount, transCount);
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

        assertEquals(4.913, allocations[0], 0.001);
        assertEquals(11.087, allocations[1], 0.001);

    }

}
