package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.neo.bind.BindConstants.REF_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.BindConstants.aminoAcidIndex;
import static com.hartwig.hmftools.neo.bind.BindCountData.INVALID_POS;
import static com.hartwig.hmftools.neo.bind.BindCountData.peptidePositionToRef;
import static com.hartwig.hmftools.neo.bind.BindCountData.refPeptidePositionToActual;

import static junit.framework.TestCase.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

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
        assertEquals(48.0, blosum.calcSequenceBlosumScore("AAA"));
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

        CalcConstants calcConstants = new CalcConstants(100, 100, 1,
                20000, 500, 500, false);

        bindCounts8.buildWeightedCounts(bindCounts, calcConstants);
        bindCounts9.buildWeightedCounts(bindCounts, calcConstants);
        bindCounts12.buildWeightedCounts(bindCounts, calcConstants);

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

        CalcConstants calcConstants = new CalcConstants(100, 100, 1,
                20000, 500, 500, false);

        BlosumMapping blosum = new BlosumMapping();

        HlaSequences hlaSequences = new HlaSequences();

        Map<String,List<String>> allelePosSequences = hlaSequences.getAllelePositionSequences();
        String seq1 = "A";
        String seq2 = "D";
        String seq3 = "E";
        allelePosSequences.put(TEST_ALLELE_1, Lists.newArrayList(seq1, seq1, seq1, seq2, seq2, seq2, seq3, seq3, seq3));
        allelePosSequences.put(TEST_ALLELE_2, Lists.newArrayList(seq1, seq2, seq3, seq3, seq2, seq2, seq3, seq3, seq3));

        bindCountsA1.buildWeightedCounts(Lists.newArrayList(bindCountsA1) , calcConstants);
        bindCountsA2.buildWeightedCounts(Lists.newArrayList(bindCountsA2), calcConstants);

        Map<String,Integer> alleleTotalCounts = Maps.newHashMap();
        alleleTotalCounts.put(TEST_ALLELE_1, pepLenA1Count);
        alleleTotalCounts.put(TEST_ALLELE_2, pepLenA2Count);

        bindCountsA1.buildFinalWeightedCounts(bindCounts, alleleTotalCounts, calcConstants, blosum, hlaSequences);
        bindCountsA2.buildFinalWeightedCounts(bindCounts, alleleTotalCounts, calcConstants, blosum, hlaSequences);

        int aaIndex = aminoAcidIndex('A');

        assertEquals(180.0, bindCountsA1.getFinalWeightedCounts()[aaIndex][0], 0.1); // same sequence
        assertEquals(101.3, bindCountsA1.getFinalWeightedCounts()[aaIndex][1], 0.1); // A & D
        assertEquals(102.5, bindCountsA1.getFinalWeightedCounts()[aaIndex][2], 0.1); // A & E
        assertEquals(105.0, bindCountsA1.getFinalWeightedCounts()[aaIndex][3], 0.1); // D & E
        assertEquals(180.0, bindCountsA1.getFinalWeightedCounts()[aaIndex][4], 0.1); // same seq

        assertEquals(450.0, bindCountsA2.getFinalWeightedCounts()[aaIndex][0], 0.1);
        assertEquals(400.2, bindCountsA2.getFinalWeightedCounts()[aaIndex][1], 0.1);
        assertEquals(406.3, bindCountsA2.getFinalWeightedCounts()[aaIndex][3], 0.1);
        assertEquals(450.0, bindCountsA2.getFinalWeightedCounts()[aaIndex][4], 0.1);
    }
}
