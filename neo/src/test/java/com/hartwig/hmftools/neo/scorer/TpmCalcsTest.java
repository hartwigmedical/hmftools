package com.hartwig.hmftools.neo.scorer;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.pointMutationInfo;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.scorer.TpmCalculator.HIGH_PROBABILITY;
import static com.hartwig.hmftools.neo.scorer.TpmCalculator.LOW_PROBABILITY;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertNotNull;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.neo.NeoEpitopeType;

import org.junit.Test;

public class TpmCalcsTest
{
    private static double[] TEST_COHORT_CANCER_TPM = { 20, 20 };
    private static double[] TEST_COHORT_PAN_CANCER_TPM = { 30, 30 };
    private static int[] TEST_PEPTIDE_LENGTH_RANGE = { 8, 8 };

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

        tpmCalc.compute("", neoDataList);

        assertEquals(7, neoData.peptides().size());

        PeptideScoreData peptideScoreData = neoData.peptides().get(0);
        assertEquals(tpms[FS_UP], peptideScoreData.expectedTpm(), 0.01);
        assertEquals(fragments, peptideScoreData.rawEffectiveTpm(), 0.01);
        assertEquals(6, peptideScoreData.effectiveTpm(), 0.01);

        // test again with 3 neoepitopes from the same variant, creating some shared and some unique peptides
        int frags1 = 2;
        double[] tpms1 = { 10, 10 };
        NeoEpitopeData neoData1 = createNeoepitope(
                pointMutationInfo("2", 1000, "A", "T"), NeoEpitopeType.MISSENSE, frags1, tpms1,
                "ABCDEFGH", "IJKL");

        int frags2 = 4;
        double[] tpms2 = { 10, 20 };
        NeoEpitopeData neoData2 = createNeoepitope(
                pointMutationInfo("2", 1000, "A", "T"), NeoEpitopeType.MISSENSE, frags2, tpms2,
                "ABCDEFGH", "IJMN");

        int frags3 = 10;
        double[] tpms3 = { 10, 50 };
        NeoEpitopeData neoData3 = createNeoepitope(
                pointMutationInfo("2", 1000, "A", "T"), NeoEpitopeType.MISSENSE, frags3, tpms3,
                "ABCDEFGH", "IJMP");

        neoDataList = Lists.newArrayList(neoData1, neoData2, neoData3);

        tpmCalc.compute("", neoDataList);

        assertEquals(4, neoData1.peptides().size());

        double tpmDownTotal = tpms1[FS_DOWN] + tpms2[FS_DOWN] + tpms3[FS_DOWN];

        PeptideScoreData peptide = neoData1.peptides().stream().filter(x -> x.Peptide.equals("BCDEFGHI")).findFirst().orElse(null);
        assertNotNull(peptide);
        assertEquals(tpms1[FS_UP], peptide.tpmUp());
        assertEquals(tpmDownTotal, peptide.tpmDownTotal());
        assertEquals(tpms1[FS_UP], peptide.expectedTpm(), 0.1);
        assertEquals(9.5, peptide.effectiveTpm(), 0.1);

        peptide = neoData1.peptides().stream().filter(x -> x.Peptide.equals("DEFGHIJK")).findFirst().orElse(null);
        assertNotNull(peptide);
        assertEquals(tpms1[FS_UP], peptide.tpmUp());
        assertEquals(tpms1[FS_DOWN], peptide.tpmDownTotal());
        double downTpm = tpms1[FS_DOWN] / tpmDownTotal;
        assertEquals(tpms1[FS_UP] * downTpm, peptide.expectedTpm(), 0.1);
        assertEquals(1.75, peptide.effectiveTpm(), 0.1);

        double tpmDownTotal2 = tpms2[FS_DOWN] + tpms3[FS_DOWN];
        peptide = neoData2.peptides().stream().filter(x -> x.Peptide.equals("DEFGHIJM")).findFirst().orElse(null);
        assertNotNull(peptide);
        assertEquals(tpms2[FS_UP], peptide.tpmUp());
        assertEquals(tpmDownTotal2, peptide.tpmDownTotal());
        assertEquals(tpms2[FS_UP] * tpmDownTotal2 / tpmDownTotal, peptide.expectedTpm(), 0.1);
        assertEquals(9.25, peptide.effectiveTpm(), 0.1);

        peptide = neoData3.peptides().stream().filter(x -> x.Peptide.equals("EFGHIJMP")).findFirst().orElse(null);
        assertNotNull(peptide);
        assertEquals(tpms3[FS_UP], peptide.tpmUp());
        assertEquals(tpms3[FS_DOWN], peptide.tpmDownTotal());
        assertEquals(tpms3[FS_UP] *  tpms3[FS_DOWN] / tpmDownTotal, peptide.expectedTpm(), 0.1);
        assertEquals(8.75, peptide.effectiveTpm(), 0.1);
    }

    private NeoEpitopeData createNeoepitope(
            final String variantInfo, final NeoEpitopeType type, int rnaFragments, final double[] tpms, final String upAAs, final String downAAs)
    {
        NeoEpitopeData neoData = new NeoEpitopeData(
                nextNeId(), type, variantInfo, GENE_ID_1, GENE_NAME_1, upAAs, "", downAAs,
                Collections.emptyList(), Collections.emptyList());

        neoData.RnaData.setCoverage(rnaFragments, rnaFragments * 2, rnaFragments * 2);
        neoData.RnaData.setExpression(tpms);
        neoData.RnaData.setCohortValues(TEST_COHORT_CANCER_TPM, TEST_COHORT_PAN_CANCER_TPM);

        return neoData;
    }

    @Test
    public void testPoissonRangeCalcs()
    {
        double reqProbLow = 0.05;
        double reqProbHigh = 0.95;

        double lowValue = TpmCalculator.calcPoissonMean(10, reqProbLow);
        double highValue = TpmCalculator.calcPoissonMean(10, reqProbHigh);

        assertEquals(17.0, lowValue, 0.1);
        assertEquals(6.3, highValue, 0.1);

        lowValue = TpmCalculator.calcPoissonMean(0, reqProbLow);
        highValue = TpmCalculator.calcPoissonMean(0, reqProbHigh);

        assertEquals(3, lowValue, 0.1);
        assertEquals(0, highValue, 0.1);

        lowValue = TpmCalculator.calcPoissonMean(1, reqProbLow);
        highValue = TpmCalculator.calcPoissonMean(1, reqProbHigh);

        assertEquals(4.7, lowValue, 0.1);
        assertEquals(0.3, highValue, 0.1);

        lowValue = TpmCalculator.calcPoissonMean(2, reqProbLow);
        highValue = TpmCalculator.calcPoissonMean(2, reqProbHigh);

        assertEquals(6.3, lowValue, 0.1);
        assertEquals(0.8, highValue, 0.1);
    }

}
