package com.hartwig.hmftools.neo.score;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.pointMutationInfo;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.neo.score.TpmCalculator.HIGH_PROBABILITY;
import static com.hartwig.hmftools.neo.score.TpmCalculator.LOW_PROBABILITY;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertNotNull;
import static junit.framework.TestCase.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.neo.NeoEpitopeType;

import org.junit.Ignore;
import org.junit.Test;

public class TpmCalcsTest
{
    private static final double[] TEST_COHORT_CANCER_TPM = { 20, 20 };
    private static final double[] TEST_COHORT_PAN_CANCER_TPM = { 30, 30 };
    private static final int[] TEST_PEPTIDE_LENGTH_RANGE = { 8, 8 };

    private static final String TRANS_NAME_1 = "TRANS_01";
    private static final String TRANS_NAME_2 = "TRANS_02";
    private static final String TRANS_NAME_3 = "TRANS_03";
    private static final String TRANS_NAME_4 = "TRANS_04";
    private static final String TRANS_NAME_5 = "TRANS_05";

    private int mNextNeId;

    public TpmCalcsTest()
    {
        mNextNeId = 0;
    }

    private int nextNeId() { return mNextNeId++; }

    @Test
    public void testEffectiveTpmCalcs()
    {
        TpmCalculator tpmCalc = new TpmCalculator(0, TEST_PEPTIDE_LENGTH_RANGE);

        List<NeoEpitopeData> neoDataList = Lists.newArrayList();

        // basic example of a missense variant creating a single neoepitope
        int fragments = 5;
        double[] tpms = { 10, 10 };

        NeoEpitopeData neoData = createNeoepitope(
                pointMutationInfo("2", 1000, "A", "T"), NeoEpitopeType.MISSENSE, fragments, tpms,
                "ABCDEFGH", "IJKLMNOP");

        neoDataList.add(neoData);

        tpmCalc.compute("", neoDataList, 2);

        assertEquals(7, neoData.peptides().size());
    }

    @Test
    public void testComplexPeptideTpmCalcs()
    {
        int peptideLenMin = 4;
        int peptideLenMax = 6;

        TpmCalculator tpmCalc = new TpmCalculator(0, new int[] { peptideLenMin, peptideLenMax });

        int fragments = 5;
        String varInfo = "1:1000-2:2000";

        double[] tpms1 = { 10, 50 };
        NeoEpitopeData neoData1 = createNeoepitope(varInfo, NeoEpitopeType.INFRAME_FUSION, fragments, tpms1, "ABC", "DEF");
        neoData1.Transcripts[FS_UP].add(TRANS_NAME_1);
        neoData1.Transcripts[FS_DOWN].add(TRANS_NAME_1);

        double[] tpms2 = { 20, 60 };
        NeoEpitopeData neoData2 = createNeoepitope(varInfo, NeoEpitopeType.INFRAME_FUSION, fragments, tpms2, "XBC", "DEF");
        neoData2.Transcripts[FS_UP].add(TRANS_NAME_2);
        neoData2.Transcripts[FS_DOWN].add(TRANS_NAME_3);

        double[] tpms3 = { 30, 70 };
        NeoEpitopeData neoData3 = createNeoepitope(varInfo, NeoEpitopeType.INFRAME_FUSION, fragments, tpms3, "ZBC", "DEG");
        neoData3.Transcripts[FS_UP].add(TRANS_NAME_3);
        neoData3.Transcripts[FS_DOWN].add(TRANS_NAME_4);

        double[] tpms4 = { 10, 80 };
        NeoEpitopeData neoData4 = createNeoepitope(varInfo, NeoEpitopeType.INFRAME_FUSION, fragments, tpms4, "ABC", "DEH");
        neoData4.Transcripts[FS_UP].add(TRANS_NAME_1); // matches 1
        neoData4.Transcripts[FS_DOWN].add(TRANS_NAME_2);

        List<NeoEpitopeData> neoDataList = Lists.newArrayList(neoData1, neoData2, neoData3, neoData4);

        tpmCalc.compute("", neoDataList, 1);

        assertEquals(6, neoData1.peptides().size());

        PeptideScoreData peptide = getPeptide("ABCD", neoDataList);
        assertEquals(10, peptide.effectiveTpm(), 0.1);

        peptide = getPeptide("ABCDEF", neoDataList);
        assertEquals(3.8, peptide.effectiveTpm(), 0.1);

        peptide = getPeptide("XBCD", neoDataList);
        assertEquals(20, peptide.effectiveTpm(), 0.1);

        peptide = getPeptide("BCDE", neoDataList);
        assertEquals(60, peptide.effectiveTpm(), 0.1);

        peptide = getPeptide("BCDEF", neoDataList);
        assertEquals(23.8, peptide.effectiveTpm(), 0.1);

        peptide = getPeptide("BCDEH", neoDataList);
        assertEquals(6.2, peptide.effectiveTpm(), 0.1);
    }

    private static PeptideScoreData getPeptide(final String peptide, final List<NeoEpitopeData> neoDataList)
    {
        for(NeoEpitopeData neoData : neoDataList)
        {
            PeptideScoreData peptideScoreData = neoData.peptides().stream().filter(x -> x.Peptide.equals(peptide)).findFirst().orElse(null);

            if(peptideScoreData != null)
                return peptideScoreData;
        }

        return null;
    }

    private NeoEpitopeData createNeoepitope(
            final String variantInfo, final NeoEpitopeType type, int rnaFragments, final double[] tpms, final String upAAs, final String downAAs)
    {
        NeoEpitopeData neoData = new NeoEpitopeData(
                nextNeId(), type, variantInfo, GENE_ID_1, GENE_NAME_1, upAAs, "", downAAs,
                Lists.newArrayList(), Lists.newArrayList(), 0, 0, 1, 1,
                0.0, 0, 0, 0, 0, 0);

        neoData.RnaData.setCoverage(rnaFragments, rnaFragments * 2, rnaFragments * 2);
        neoData.RnaData.setExpression(tpms);
        neoData.RnaData.setCohortValues(TEST_COHORT_CANCER_TPM, TEST_COHORT_PAN_CANCER_TPM);

        return neoData;
    }
    
    @Test
    public void testPoissonRangeCalcs()
    {
        double reqProbLow = LOW_PROBABILITY;
        double reqProbHigh = HIGH_PROBABILITY;

        double lowValue = TpmCalculator.calcPoissonMean(10, reqProbLow);
        double highValue = TpmCalculator.calcPoissonMean(10, reqProbHigh);

        assertEquals(13.0, lowValue, 0.1);
        assertEquals(8.6, highValue, 0.1);

        lowValue = TpmCalculator.calcPoissonMean(0, reqProbLow);
        highValue = TpmCalculator.calcPoissonMean(0, reqProbHigh);

        assertEquals(1.39, lowValue, 0.1);
        assertEquals(0, highValue, 0.1);

        lowValue = TpmCalculator.calcPoissonMean(1, reqProbLow);
        highValue = TpmCalculator.calcPoissonMean(1, reqProbHigh);

        assertEquals(2.7, lowValue, 0.1);
        assertEquals(0.876, highValue, 0.1);

        lowValue = TpmCalculator.calcPoissonMean(2, reqProbLow);
        highValue = TpmCalculator.calcPoissonMean(2, reqProbHigh);

        assertEquals(3.9, lowValue, 0.1);
        assertEquals(1.7, highValue, 0.1);
    }

}
