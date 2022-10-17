package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.MatrixUtils.clear;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_ALLELE;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PEPTIDE_LEN;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_ALLELE_MOTIF_WEIGHT_FACTOR;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_ALLELE_MOTIF_WEIGHT_MAX;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_PEP_LEN_WEIGHT_FACTOR;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_PEP_LEN_WEIGHT_MAX;
import static com.hartwig.hmftools.neo.bind.BindConstants.aminoAcidIndex;
import static com.hartwig.hmftools.neo.bind.NoiseModel.calcExpected;
import static com.hartwig.hmftools.neo.bind.PosWeightModel.INVALID_POS;
import static com.hartwig.hmftools.neo.bind.PosWeightModel.peptidePosToRef;
import static com.hartwig.hmftools.neo.bind.PosWeightModel.refPeptidePosToActual;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.VectorUtils;
import com.hartwig.hmftools.neo.utils.AminoAcidFrequency;

import org.junit.Test;

public class BindModelTest
{
    private static final String TEST_ALLELE_1 = "A0101";
    private static final String TEST_ALLELE_2 = "A0201";
    private static final String TEST_ALLELE_3 = "A0301";
    private static final String TEST_ALLELE_4 = "A0401";

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
        assertEquals(0, peptidePosToRef(8, 0));
        assertEquals(4, peptidePosToRef(8, 4));
        assertEquals(9, peptidePosToRef(8, 5));
        assertEquals(11, peptidePosToRef(8, 7));

        assertEquals(8, peptidePosToRef(9, 5));
        assertEquals(11, peptidePosToRef(9, 8));
        assertEquals(11, peptidePosToRef(10, 9));
        assertEquals(11, peptidePosToRef(11, 10));
        assertEquals(11, peptidePosToRef(12, 11));

        assertEquals(5, refPeptidePosToActual(8, 9));
        assertEquals(7, refPeptidePosToActual(8, 11));
        assertEquals(11, refPeptidePosToActual(12, 11));

        assertEquals(INVALID_POS, refPeptidePosToActual(8, 8));
    }


    @Test
    public void testPeptideLengthWeights()
    {
        BindCountData bindCounts8 = new BindCountData(TEST_ALLELE_1, 8);
        BindCountData bindCounts9 = new BindCountData(TEST_ALLELE_1, 9);
        BindCountData bindCounts10 = new BindCountData(TEST_ALLELE_1, 10);
        BindCountData bindCounts11 = new BindCountData(TEST_ALLELE_1, 11);

        int testPos = 0;
        int aa = 0;
        bindCounts8.getBindCounts()[aa][testPos] = 30;
        bindCounts8.getBindCounts()[aa+1][testPos] = 300 - bindCounts8.getBindCounts()[aa][testPos];

        bindCounts9.getBindCounts()[aa][testPos] = 80;
        bindCounts9.getBindCounts()[aa+1][testPos] = 800 - bindCounts9.getBindCounts()[aa][testPos];

        bindCounts10.getBindCounts()[aa][testPos] = 500;
        bindCounts10.getBindCounts()[aa+1][testPos] = 1000 - bindCounts10.getBindCounts()[aa][testPos];

        bindCounts11.getBindCounts()[aa][testPos] = 5;
        bindCounts11.getBindCounts()[aa+1][testPos] = 50 - bindCounts11.getBindCounts()[aa][testPos];

        List<BindCountData> bindCounts = Lists.newArrayList(bindCounts8, bindCounts9, bindCounts10, bindCounts11);

        CalcConstants calcConstants = new CalcConstants(0.1, 200,
                0, 0, 0, 0);

        PosWeightModel pwModel = new PosWeightModel(calcConstants, new HlaSequences());

        pwModel.buildWeightedCounts(bindCounts8, bindCounts);
        pwModel.buildWeightedCounts(bindCounts9, bindCounts);
        pwModel.buildWeightedCounts(bindCounts10, bindCounts);
        pwModel.buildWeightedCounts(bindCounts11, bindCounts);

        assertEquals(35.9, bindCounts8.getWeightedCounts()[aa][testPos], 0.1);
        assertEquals(85.4, bindCounts9.getWeightedCounts()[aa][testPos], 0.1);
        assertEquals(501.2, bindCounts10.getWeightedCounts()[aa][testPos], 0.1);
        assertEquals(10.8, bindCounts11.getWeightedCounts()[aa][testPos], 0.1);

        // check position translation
        clearAll(bindCounts8);
        clearAll(bindCounts11);

        bindCounts8.getBindCounts()[aa][5] = 100; // 5 -> 8 in 11

        bindCounts11.getBindCounts()[aa][testPos] = 0;
        bindCounts11.getBindCounts()[aa+1][testPos] = 0;

        bindCounts11.getBindCounts()[aa][6] = 100; // no impact
        bindCounts11.getBindCounts()[aa][7] = 100; // no impact
        bindCounts11.getBindCounts()[aa][8] = 100;

        bindCounts = Lists.newArrayList(bindCounts8, bindCounts11);

        pwModel.buildWeightedCounts(bindCounts8, bindCounts);
        pwModel.buildWeightedCounts(bindCounts11, bindCounts);

        assertEquals(101.0, bindCounts8.getWeightedCounts()[aa][5], 0.1);
        assertEquals(101.0, bindCounts11.getWeightedCounts()[aa][8], 0.1);

        assertEquals(100.0, bindCounts11.getWeightedCounts()[aa][6], 0.1);
        assertEquals(100.0, bindCounts11.getWeightedCounts()[aa][7], 0.1);
    }

    @Test
    public void testAlleleMotifWeights()
    {
        // a test to match and demonstrate the documentation example calcs
        List<BindCountData> allBindCounts = loadTestBindCounts("/allele_pos_counts_1.csv", 0);

        Set<String> alleles = Sets.newHashSet();
        allBindCounts.forEach(x -> alleles.add(x.Allele));
        assertEquals(4, alleles.size());

        HlaSequences hlaSequences = new HlaSequences();

        Map<String,List<String>> allelePosSequences = hlaSequences.getAllelePositionSequences();
        String seqOther = "AAA";
        List<String> otherSequences = Lists.newArrayList(seqOther, seqOther, seqOther, seqOther, seqOther, seqOther, seqOther, seqOther);

        List<String> sequences = Lists.newArrayList("AHV");
        sequences.addAll(otherSequences);
        allelePosSequences.put(TEST_ALLELE_1, sequences);

        sequences = Lists.newArrayList("AHA");
        sequences.addAll(otherSequences);
        allelePosSequences.put(TEST_ALLELE_2, sequences);

        sequences = Lists.newArrayList("AKA");
        sequences.addAll(otherSequences);
        allelePosSequences.put(TEST_ALLELE_3, sequences);

        sequences = Lists.newArrayList("AHA");
        sequences.addAll(otherSequences);
        allelePosSequences.put(TEST_ALLELE_4, sequences);

        CalcConstants calcConstants = new CalcConstants(0.2, 200,
                0.2, 100, 0, 0);

        PosWeightModel pwModel = new PosWeightModel(calcConstants, hlaSequences);

        for(String allele : alleles)
        {
            List<BindCountData> alleleBindCounts = allBindCounts.stream().filter(x -> x.Allele.equals(allele)).collect(Collectors.toList());
            alleleBindCounts.forEach(x -> pwModel.buildWeightedCounts(x, alleleBindCounts));
            pwModel.buildPositionAdjustedTotals(alleleBindCounts);
        }

        BindCountData bindCounts9 = allBindCounts.stream()
                .filter(x -> x.PeptideLength == 9 && x.Allele.equals(TEST_ALLELE_4)).findFirst().orElse(null);
        assertTrue(bindCounts9 != null);

        int testPos = 0;
        int testAa = 0;
        assertEquals(854.0, bindCounts9.getAdjustedPosTotals()[testPos], 0.1);

        // now test allele-motif weighting
        for(String allele : alleles)
        {
            List<BindCountData> alleleBindCounts = allBindCounts.stream().filter(x -> x.Allele.equals(allele)).collect(Collectors.toList());

            for(BindCountData bindCounts : alleleBindCounts)
            {
                List<BindCountData> pepLenBindCounts =
                        allBindCounts.stream().filter(x -> x.PeptideLength == bindCounts.PeptideLength).collect(Collectors.toList());

                pwModel.buildFinalWeightedCounts(bindCounts, pepLenBindCounts);
            }
        }

        assertEquals(97.9, bindCounts9.getFinalWeightedCounts()[testAa][testPos], 0.1);

        BindCountData bindCounts8 = allBindCounts.stream()
                .filter(x -> x.PeptideLength == 8 && x.Allele.equals(TEST_ALLELE_4)).findFirst().orElse(null);

        assertEquals(46.0, bindCounts8.getFinalWeightedCounts()[testAa][testPos], 0.1);

        BindCountData bindCounts10 = allBindCounts.stream()
                .filter(x -> x.PeptideLength == 10 && x.Allele.equals(TEST_ALLELE_4)).findFirst().orElse(null);

        assertEquals(505.1, bindCounts10.getFinalWeightedCounts()[testAa][testPos], 0.1);

        BindCountData bindCounts11 = allBindCounts.stream()
                .filter(x -> x.PeptideLength == 11 && x.Allele.equals(TEST_ALLELE_4)).findFirst().orElse(null);

        assertEquals(18.4, bindCounts11.getFinalWeightedCounts()[testAa][testPos], 0.1);
    }

    private static void clearAll(final BindCountData bindCounts)
    {
        clear(bindCounts.getBindCounts());
        clear(bindCounts.getWeightedCounts());
        clear(bindCounts.getFinalWeightedCounts());
        VectorUtils.clear(bindCounts.getAdjustedPosTotals());
    }

    private static List<BindCountData> loadTestBindCounts(final String bindCountsFile, int targetPosIndex)
    {
        List<BindCountData> bindCountDataList = Lists.newArrayList();

        final List<String> lines = new BufferedReader(new InputStreamReader(
                BindModelTest.class.getResourceAsStream(bindCountsFile))).lines().collect(Collectors.toList());

        final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);
        lines.remove(0);

        int alleleIndex = fieldsIndexMap.get(FLD_ALLELE);
        int pepLenIndex = fieldsIndexMap.get(FLD_PEPTIDE_LEN);
        int aaCountIndex = fieldsIndexMap.get("AminoAcidCount");
        int posTotalIndex = fieldsIndexMap.get("PositionTotal");

        for(String line : lines)
        {
            final String[] values = line.split(DELIMITER, -1);

            BindCountData bindCounts = new BindCountData(values[alleleIndex], Integer.parseInt(values[pepLenIndex]));
            int aaCount = Integer.parseInt(values[aaCountIndex]);
            int posTotal = Integer.parseInt(values[posTotalIndex]);
            bindCounts.getBindCounts()[0][targetPosIndex] = aaCount;
            bindCounts.getBindCounts()[1][targetPosIndex] = posTotal - aaCount;
            bindCountDataList.add(bindCounts);
        }

        return bindCountDataList;
    }

    /*
    @Test
    public void testAlleleMotifWeights()
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
    */

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
        return new CalcConstants(DEFAULT_PEP_LEN_WEIGHT_FACTOR, DEFAULT_PEP_LEN_WEIGHT_MAX,
                DEFAULT_ALLELE_MOTIF_WEIGHT_FACTOR, DEFAULT_ALLELE_MOTIF_WEIGHT_MAX, 0, 0);
    }
}
