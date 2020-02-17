package com.hartwig.hmftools.svtools.rna_expression;

import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.ExpectedExpressionRates.FL_SIZE;
import static com.hartwig.hmftools.svtools.rna_expression.ExpectedExpressionRates.UNSPLICED_ID;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TC_LONG;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TC_SHORT;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TC_SPLICED;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TC_UNSPLICED;
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

import org.junit.Test;

public class RnaExpectedRates
{
    @Test
    public void testRegionMatching()
    {
        RnaExpConfig config = new RnaExpConfig();
        int fragmentLength = 100;
        config.ReadLength = 20;

        ExpectedExpressionRates eeRates = new ExpectedExpressionRates(config);

        eeRates.setFragmentLengthData(fragmentLength, 1);

        TranscriptData transData = new TranscriptData(1, "TRANS01", "GENE01", true, (byte)1,
                0, 1000, null,null, "");

        transData.exons().add(new ExonData(1, 100, 200, 1, -1, -1));
        transData.exons().add(new ExonData(1, 300, 400, 2, -1, -1));
        transData.exons().add(new ExonData(1, 440, 449, 3, -1, -1));
        transData.exons().add(new ExonData(1, 460, 469, 4, -1, -1));
        transData.exons().add(new ExonData(1, 600, 800, 5, -1, -1));

        // fully contained fragment
        long startPos = 100;
        List<long[]> readRegions = Lists.newArrayList();
        List<long[]> spliceJunctions = Lists.newArrayList();
        int matchType = eeRates.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertEquals(TC_SHORT, matchType);
        assertEquals(2, readRegions.size());
        assertEquals(100, readRegions.get(0)[SE_START]);
        assertEquals(119, readRegions.get(0)[SE_END]);
        assertEquals(180, readRegions.get(1)[SE_START]);
        assertEquals(199, readRegions.get(1)[SE_END]);
        assertTrue(eeRates.readsSupportFragment(transData, readRegions, matchType, spliceJunctions));

        // test for another transcript
        TranscriptData transData2 = new TranscriptData(2, "TRANS02", "GENE01", true, (byte)1,
                0, 1000, null,null, "");

        transData2.exons().add(new ExonData(2, 90, 210, 1, -1, -1));
        transData2.exons().add(new ExonData(2, 300, 400, 2, -1, -1));
        assertTrue(eeRates.readsSupportFragment(transData2, readRegions, matchType, spliceJunctions));

        transData2.exons().clear();
        transData2.exons().add(new ExonData(2, 150, 250, 1, -1, -1));
        transData2.exons().add(new ExonData(2, 300, 400, 2, -1, -1));
        assertFalse(eeRates.readsSupportFragment(transData2, readRegions, matchType, spliceJunctions));

        // 2 fully contained reads but in 2 exons
        startPos = 371;
        matchType = eeRates.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertEquals(TC_LONG, matchType);
        assertEquals(2, readRegions.size());
        assertEquals(371, readRegions.get(0)[SE_START]);
        assertEquals(390, readRegions.get(0)[SE_END]);
        assertEquals(631, readRegions.get(1)[SE_START]);
        assertEquals(650, readRegions.get(1)[SE_END]);
        assertTrue(eeRates.readsSupportFragment(transData, readRegions, matchType, spliceJunctions));

        transData2.exons().clear();
        transData2.exons().add(new ExonData(2, 250, 400, 1, -1, -1));
        transData2.exons().add(new ExonData(2, 450, 550, 2, -1, -1));
        transData2.exons().add(new ExonData(2, 600, 700, 3, -1, -1));
        assertTrue(eeRates.readsSupportFragment(transData2, readRegions, matchType, spliceJunctions));

        transData2.exons().clear();
        transData2.exons().add(new ExonData(2, 250, 300, 1, -1, -1));
        transData2.exons().add(new ExonData(2, 350, 700, 2, -1, -1));
        assertTrue(eeRates.readsSupportFragment(transData2, readRegions, matchType, spliceJunctions));

        // within an exon then spanning a junction
        startPos = 150;
        fragmentLength = 32 + 85 + 2 * config.ReadLength;
        eeRates.setFragmentLengthData(fragmentLength, 1);

        matchType = eeRates.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertEquals(TC_SPLICED, matchType);
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

        assertEquals(TC_SPLICED, matchType);
        assertEquals(3, readRegions.size());
        assertEquals(191, readRegions.get(0)[SE_START]);
        assertEquals(200, readRegions.get(0)[SE_END]);
        assertEquals(300, readRegions.get(1)[SE_START]);
        assertEquals(309, readRegions.get(1)[SE_END]);
        assertEquals(350, readRegions.get(2)[SE_START]);
        assertEquals(369, readRegions.get(2)[SE_END]);
        assertTrue(eeRates.readsSupportFragment(transData, readRegions, matchType, spliceJunctions));

        transData2.exons().clear();
        transData2.exons().add(new ExonData(2, 150, 200, 1, -1, -1));
        transData2.exons().add(new ExonData(2, 300, 320, 2, -1, -1));
        transData2.exons().add(new ExonData(2, 340, 380, 3, -1, -1));
        assertTrue(eeRates.readsSupportFragment(transData2, readRegions, matchType, spliceJunctions));

        // invalid since an exon is skipped in the splicing read
        transData2.exons().clear();
        transData2.exons().add(new ExonData(2, 150, 200, 1, -1, -1));
        transData2.exons().add(new ExonData(2, 210, 290, 2, -1, -1));
        transData2.exons().add(new ExonData(2, 300, 450, 3, -1, -1));
        assertFalse(eeRates.readsSupportFragment(transData2, readRegions, matchType, spliceJunctions));

        // 2 sets of exon junctions
        startPos = 191;
        fragmentLength = 81 + 5 + 2 * config.ReadLength;
        eeRates.setFragmentLengthData(fragmentLength, 1);

        matchType = eeRates.generateImpliedFragment(transData, startPos, readRegions, spliceJunctions);

        assertEquals(TC_SPLICED, matchType);
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
    }

    @Test
    public void testSingleTranscriptCounts()
    {
        RnaExpConfig config = new RnaExpConfig();
        config.ReadLength = 10;
        config.ExpRateFragmentLengths.add(new int[] {30, 1});

        ExpectedExpressionRates eeRates = new ExpectedExpressionRates(config);

        String geneId = "GENE01";

        EnsemblGeneData geneData = new EnsemblGeneData(geneId, geneId, "1", (byte)1, 100, 1000, "");

        String transId = "TRANS01";

        TranscriptData transData = new TranscriptData(1, transId, geneId, true, (byte)1,
                100, 414, null,null, "");

        transData.exons().add(new ExonData(1, 100, 158, 1, -1, -1));
        transData.exons().add(new ExonData(1, 228, 286, 2, -1, -1));
        transData.exons().add(new ExonData(1, 356, 414, 3, -1, -1));

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

        List<String> tranIds = Lists.newArrayList(transId);

        TranscriptComboData tcData = findMatchingData(tranIds, tcDataList);
        assertTrue(tcData != null);
        assertEquals(90, tcData.getCount(TC_SHORT));
        assertEquals(58, tcData.getCount(TC_SPLICED));

        tcDataList = transComboData.get(UNSPLICED_ID);
        assertEquals(2, tcDataList.size());

        tcData = findMatchingData(tranIds, tcDataList);
        assertTrue(tcData != null);
        assertEquals(87, tcData.getCount(TC_SHORT));
        assertEquals(0, tcData.getCount(TC_SPLICED));

        tranIds.clear();
        tcData = findMatchingData(tranIds, tcDataList);
        assertTrue(tcData != null);
        assertEquals(196, tcData.getCount(TC_UNSPLICED));
        assertEquals(tcData.getCount(TC_UNSPLICED), tcData.totalCount());
    }

    @Test
    public void testMultipleTranscriptCounts()
    {
        RnaExpConfig config = new RnaExpConfig();
        config.ExpRateFragmentLengths.add(new int[] {30, 1});
        config.ReadLength = 10;

        ExpectedExpressionRates eeRates = new ExpectedExpressionRates(config);

        String geneId = "GENE01";

        EnsemblGeneData geneData = new EnsemblGeneData(geneId, geneId, "1", (byte)1, 100, 1000, "");

        String transId1 = "TRANS01";

        TranscriptData transData1 = new TranscriptData(1, transId1, geneId, true, (byte)1,
                100, 600, null,null, "");

        transData1.exons().add(new ExonData(1, 100, 200, 1, -1, -1));
        transData1.exons().add(new ExonData(1, 300, 400, 2, -1, -1));
        transData1.exons().add(new ExonData(1, 500, 600, 3, -1, -1));

        String transId2 = "TRANS02";

        TranscriptData transData2 = new TranscriptData(1, transId2, geneId, true, (byte)1,
                150, 1000, null,null, "");

        transData2.exons().add(new ExonData(1, 150, 200, 1, -1, -1));
        transData2.exons().add(new ExonData(1, 300, 450, 2, -1, -1));
        transData2.exons().add(new ExonData(1, 700, 800, 3, -1, -1));
        transData2.exons().add(new ExonData(1, 900, 1000, 4, -1, -1));

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

        List<String> tranIds = Lists.newArrayList(transId1);

        TranscriptComboData tcData = findMatchingData(tranIds, tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.getCount(TC_SHORT) > 0);
        assertTrue(tcData.getCount(TC_SPLICED) > 0);
        assertEquals(tcData.getCount(TC_SPLICED) + tcData.getCount(TC_SHORT), tcData.totalCount());

        tranIds = Lists.newArrayList(transId1, transId2);
        tcData = findMatchingData(tranIds, tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.getCount(TC_SHORT) > 0);
        assertTrue(tcData.getCount(TC_SPLICED) > 0);

        tcDataList = transComboData.get(transId2);
        assertEquals(2, tcDataList.size());

        tranIds = Lists.newArrayList(transId2);

        tcData = findMatchingData(tranIds, tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.getCount(TC_SHORT) > 0);
        assertTrue(tcData.getCount(TC_SPLICED) > 0);

        tranIds = Lists.newArrayList(transId1, transId2);
        tcData = findMatchingData(tranIds, tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.getCount(TC_SHORT) > 0);
        assertTrue(tcData.getCount(TC_SPLICED) > 0);

        tcDataList = transComboData.get(UNSPLICED_ID);
        assertEquals(4, tcDataList.size());

        tranIds.clear();
        tcData = findMatchingData(tranIds, tcDataList);
        assertTrue(tcData != null);
        assertTrue(tcData.getCount(TC_UNSPLICED) > 0);
        assertEquals(tcData.getCount(TC_UNSPLICED), tcData.totalCount());
    }

}
