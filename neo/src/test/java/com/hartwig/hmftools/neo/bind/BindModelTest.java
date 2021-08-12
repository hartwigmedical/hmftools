package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.neo.bind.BindConstants.REF_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.BindConstants.aminoAcidIndex;
import static com.hartwig.hmftools.neo.bind.NoiseModel.calcExpected;
import static com.hartwig.hmftools.neo.bind.PosWeightModel.INVALID_POS;
import static com.hartwig.hmftools.neo.bind.PosWeightModel.peptidePositionToRef;
import static com.hartwig.hmftools.neo.bind.PosWeightModel.refPeptidePositionToActual;

import static junit.framework.TestCase.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.neo.utils.AminoAcidFrequency;

import org.junit.Test;

public class BindModelTest
{
    private static final String TEST_ALLELE_1 = "A0101";
    private static final String TEST_ALLELE_2 = "A0102";

    @Test
    public void testBindCounts()
    {
        BindCountData bindCounts = new BindCountData(TEST_ALLELE_1, 8);

        char AA_1 = 'E';
        char AA_2 = 'D';
        String peptide1 = "ACD" + AA_1 + "FGHI";
        String peptide2 = "ACD" + AA_2 + "FGHI";

        for(int i = 0; i < 20; ++i)
        {
            bindCounts.processBindingPeptide(peptide1);
        }

        for(int i = 0; i < 50; ++i)
        {
            bindCounts.processBindingPeptide(peptide2);
        }

        int aa1Index = aminoAcidIndex(AA_1);
        int aa2Index = aminoAcidIndex(AA_2);

        assertEquals(70, bindCounts.totalBindCount());

        int pos = 3;
        assertEquals(20.0, bindCounts.getBindCounts()[aa1Index][pos]);
        assertEquals(50.0, bindCounts.getBindCounts()[aa2Index][pos]);
    }

    @Test
    public void testPeptideLengthConversion()
    {
        assertEquals(0, peptidePositionToRef(REF_PEPTIDE_LENGTH, 8, 0));
        assertEquals(4, peptidePositionToRef(REF_PEPTIDE_LENGTH, 8, 4));
        assertEquals(9, peptidePositionToRef(REF_PEPTIDE_LENGTH, 8, 5));
        assertEquals(11, peptidePositionToRef(REF_PEPTIDE_LENGTH, 8, 7));

        assertEquals(8, peptidePositionToRef(REF_PEPTIDE_LENGTH, 9, 5));
        assertEquals(11, peptidePositionToRef(REF_PEPTIDE_LENGTH, 9, 8));
        assertEquals(11, peptidePositionToRef(REF_PEPTIDE_LENGTH, 10, 9));
        assertEquals(11, peptidePositionToRef(REF_PEPTIDE_LENGTH, 11, 10));
        assertEquals(11, peptidePositionToRef(REF_PEPTIDE_LENGTH, 12, 11));

        assertEquals(5, refPeptidePositionToActual(REF_PEPTIDE_LENGTH, 8, 9));
        assertEquals(7, refPeptidePositionToActual(REF_PEPTIDE_LENGTH, 8, 11));
        assertEquals(11, refPeptidePositionToActual(REF_PEPTIDE_LENGTH, 12, 11));

        assertEquals(INVALID_POS, refPeptidePositionToActual(REF_PEPTIDE_LENGTH, 8, 8));
    }

    @Test
    public void testBlosumCalcs()
    {
        BlosumMapping blosum = new BlosumMapping();

        assertEquals(16.0, blosum.calcSequenceBlosumScore("A"));
        assertEquals(4096.0, blosum.calcSequenceBlosumScore("AAA"));
        assertEquals(128.0, blosum.calcSequenceBlosumScore("Y"));

        assertEquals(4.0, blosum.calcSequenceBlosumScore("D", "E"));
        assertEquals(4.0, blosum.calcSequenceBlosumScore("E", "D"));

        assertEquals(0.125, blosum.calcSequenceBlosumScore("W", "A"), 0.001);
        assertEquals(0.125, blosum.calcSequenceBlosumScore("A", "W"), 0.001);
    }

    @Test
    public void testPeptideLengthWeights()
    {
        BindCountData bindCounts8 = new BindCountData(TEST_ALLELE_1, 8);

        int pos = 3;
        String peptide8 = "AAAAAAAA";

        int pepLen8Count = 100;
        for(int i = 0; i < pepLen8Count; ++i)
        {
            bindCounts8.processBindingPeptide(peptide8);
        }

        int aaIndex = aminoAcidIndex('A');

        assertEquals(pepLen8Count, bindCounts8.totalBindCount());
        assertEquals(100.0, bindCounts8.getBindCounts()[aaIndex][pos]);

        BindCountData bindCounts9 = new BindCountData(TEST_ALLELE_1, 9);

        String peptide9 = "AAAAAAAAA";

        int pepLen9Count = 400;
        for(int i = 0; i < pepLen9Count; ++i)
        {
            bindCounts9.processBindingPeptide(peptide9);
        }

        assertEquals(400.0, bindCounts9.getBindCounts()[aaIndex][pos]);

        BindCountData bindCounts12 = new BindCountData(TEST_ALLELE_1, 12);

        String peptide12 = "AAAAAAAAAAAA";

        int pepLen12Count = 200;
        for(int i = 0; i < pepLen12Count; ++i)
        {
            bindCounts12.processBindingPeptide(peptide12);
        }

        List<BindCountData> bindCounts = Lists.newArrayList(bindCounts8, bindCounts9, bindCounts12);

        CalcConstants calcConstants = getTestCalcConstants();

        PosWeightModel pwModel = new PosWeightModel(calcConstants, new HlaSequences());

        pwModel.buildWeightedCounts(bindCounts8, bindCounts);
        pwModel.buildWeightedCounts(bindCounts9, bindCounts);
        pwModel.buildWeightedCounts(bindCounts12, bindCounts);

        assertEquals(325, bindCounts8.getWeightedCounts()[aaIndex][0], 0.1);
        assertEquals(325, bindCounts8.getWeightedCounts()[aaIndex][7], 0.1);

        assertEquals(433.3, bindCounts9.getWeightedCounts()[aaIndex][0], 0.1);
        assertEquals(433.3, bindCounts9.getWeightedCounts()[aaIndex][4], 0.1);
        assertEquals(413.3, bindCounts9.getWeightedCounts()[aaIndex][5], 0.1);
        assertEquals(433.3, bindCounts9.getWeightedCounts()[aaIndex][8], 0.1);

        assertEquals(252.8, bindCounts12.getWeightedCounts()[aaIndex][0], 0.1);
        assertEquals(252.8, bindCounts12.getWeightedCounts()[aaIndex][4], 0.1);
        assertEquals(200.0, bindCounts12.getWeightedCounts()[aaIndex][5]);
        assertEquals(200.0, bindCounts12.getWeightedCounts()[aaIndex][7]);
        assertEquals(244.4, bindCounts12.getWeightedCounts()[aaIndex][8], 0.1);
        assertEquals(252.8, bindCounts12.getWeightedCounts()[aaIndex][11], 0.1);
    }

    @Test
    public void setTestAlleleMotifWeights()
    {
        BindCountData bindCountsA1 = new BindCountData(TEST_ALLELE_1, 9);

        String peptide = "AAAAAAAAA";

        int pepLenA1Count = 100;
        for(int i = 0; i < pepLenA1Count; ++i)
        {
            bindCountsA1.processBindingPeptide(peptide);
        }

        BindCountData bindCountsA2 = new BindCountData(TEST_ALLELE_2, 9);

        int pepLenA2Count = 400;
        for(int i = 0; i < pepLenA2Count; ++i)
        {
            bindCountsA2.processBindingPeptide(peptide);
        }

        List<BindCountData> bindCounts = Lists.newArrayList(bindCountsA1, bindCountsA2);

        CalcConstants calcConstants = getTestCalcConstants();

        HlaSequences hlaSequences = new HlaSequences();

        Map<String,List<String>> allelePosSequences = hlaSequences.getAllelePositionSequences();
        String seq1 = "A";
        String seq2 = "D";
        String seq3 = "E";
        allelePosSequences.put(TEST_ALLELE_1, Lists.newArrayList(seq1, seq1, seq1, seq2, seq2, seq2, seq3, seq3, seq3));
        allelePosSequences.put(TEST_ALLELE_2, Lists.newArrayList(seq1, seq2, seq3, seq3, seq2, seq2, seq3, seq3, seq3));

        PosWeightModel pwModel = new PosWeightModel(calcConstants, hlaSequences);

        pwModel.buildWeightedCounts(bindCountsA1, Lists.newArrayList(bindCountsA1));
        pwModel.buildWeightedCounts(bindCountsA2, Lists.newArrayList(bindCountsA2));

        Map<String,Integer> alleleTotalCounts = Maps.newHashMap();
        alleleTotalCounts.put(TEST_ALLELE_1, pepLenA1Count);
        alleleTotalCounts.put(TEST_ALLELE_2, pepLenA2Count);

        pwModel.buildFinalWeightedCounts(bindCountsA1, bindCounts, alleleTotalCounts);
        pwModel.buildFinalWeightedCounts(bindCountsA2, bindCounts, alleleTotalCounts);

        int aaIndex = aminoAcidIndex('A');

        assertEquals(300, bindCountsA1.getFinalWeightedCounts()[aaIndex][0], 0.1); // same sequence
        assertEquals(103.1, bindCountsA1.getFinalWeightedCounts()[aaIndex][1], 0.1); // A & D
        assertEquals(106.3, bindCountsA1.getFinalWeightedCounts()[aaIndex][2], 0.1); // A & E
        assertEquals(112.5, bindCountsA1.getFinalWeightedCounts()[aaIndex][3], 0.1); // D & E
        assertEquals(300.0, bindCountsA1.getFinalWeightedCounts()[aaIndex][4], 0.1); // same seq

        assertEquals(420.0, bindCountsA2.getFinalWeightedCounts()[aaIndex][0], 0.1);
        assertEquals(400.1, bindCountsA2.getFinalWeightedCounts()[aaIndex][1], 0.1);
        assertEquals(402.5, bindCountsA2.getFinalWeightedCounts()[aaIndex][3], 0.1);
        assertEquals(420.0, bindCountsA2.getFinalWeightedCounts()[aaIndex][4], 0.1);
    }

    @Test
    public void testNoiseModel()
    {
        // first the raw probability calcs
        int totalBinds = 100;
        double reqProbability = 0.05;
        double expected = calcExpected(totalBinds, 1, 10, reqProbability);
        assertEquals(4.7, expected, 0.1);

        expected = calcExpected(totalBinds, 0, 10, reqProbability);
        assertEquals(3.0, expected, 0.1);

        expected = calcExpected(totalBinds, 5, 15, reqProbability);
        assertEquals(10.3, expected, 0.1);

        // now test using cache and amino-acid frequency

        // initially with a zero weight to fully use expected
        NoiseModel noiseModel = new NoiseModel(new AminoAcidFrequency(), 0.05, 0.0);

        // L is about 10%
        expected = noiseModel.getExpected('L', 100, 1);
        assertEquals(4.7, expected, 0.1);

        // value cached for a similar request
        expected = noiseModel.getExpected('L', 100, 1);
        assertEquals(1, noiseModel.cacheSize());

        // test larger values with rounding
        expected = noiseModel.getExpected('L', 1504, 47);
        assertEquals(60, expected, 1.0);

        noiseModel.getExpected('L', 1504, 51);
        assertEquals(3, noiseModel.cacheSize());

        noiseModel.getExpected('L', 1506, 51); // observed is rounded to as previous call
        assertEquals(3, noiseModel.cacheSize());

        // and with roundng the observed
        expected = noiseModel.getExpected('L', 10100, 802); // observed is rounded to as previous call
        assertEquals(847, expected, 1.0);
        assertEquals(4, noiseModel.cacheSize());

        noiseModel.getExpected('L', 10099, 799); // observed is rounded to as previous call
        noiseModel.getExpected('L', 10103, 803); // observed is rounded to as previous call
        assertEquals(4, noiseModel.cacheSize());

        // test weight
        noiseModel = new NoiseModel(new AminoAcidFrequency(), 0.05, 0.5);
        expected = noiseModel.getExpected('L', 10100, 802); // observed is rounded to as previous call
        assertEquals(824, expected, 1.0);
    }

    private static CalcConstants getTestCalcConstants()
    {
        return new CalcConstants(100, 100, 1, 0, 0, 0);
    }
}
