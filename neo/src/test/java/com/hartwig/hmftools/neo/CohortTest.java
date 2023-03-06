package com.hartwig.hmftools.neo;

import static com.hartwig.hmftools.neo.scorer.AlleleCoverage.EXPECTED_ALLELE_COUNT;

import static junit.framework.TestCase.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.neo.scorer.PeptidePredictionData;
import com.hartwig.hmftools.neo.scorer.TpmCalculator;

import org.junit.Test;

public class CohortTest
{
    private static PeptidePredictionData make(String allele, int neId, String peptide, double affinity, double presentation)
    {
        PeptidePredictionData predData = new PeptidePredictionData(neId, allele, peptide);
        return predData;
    }

    @Test
    public void testPredictionAlleles()
    {
        List<PeptidePredictionData> predictions = Lists.newArrayList();

        // test all homozygous
        predictions.add(make("A0101", 1, "ABC", 1, 0.95));
        predictions.add(make("B0101", 1, "ABC", 1, 0.95));
        predictions.add(make("C0101", 1, "ABC", 1, 0.95));

        PeptidePredictionData.expandHomozygous(predictions);
        assertEquals(EXPECTED_ALLELE_COUNT, predictions.size());
        assertEquals(2, predictions.stream().filter(x -> x.Allele.equals("A0101")).count());
        assertEquals(2, predictions.stream().filter(x -> x.Allele.equals("B0101")).count());
        assertEquals(2, predictions.stream().filter(x -> x.Allele.equals("C0101")).count());

        // test each of only 1 gene being homozygous
        predictions.remove(0);
        PeptidePredictionData.expandHomozygous(predictions);
        assertEquals(EXPECTED_ALLELE_COUNT, predictions.size());
        assertEquals(2, predictions.stream().filter(x -> x.Allele.equals("A0101")).count());

        predictions.remove(2);
        PeptidePredictionData.expandHomozygous(predictions);
        assertEquals(EXPECTED_ALLELE_COUNT, predictions.size());
        assertEquals(2, predictions.stream().filter(x -> x.Allele.equals("B0101")).count());

        predictions.remove(4);
        PeptidePredictionData.expandHomozygous(predictions);
        assertEquals(EXPECTED_ALLELE_COUNT, predictions.size());
        assertEquals(2, predictions.stream().filter(x -> x.Allele.equals("C0101")).count());
    }

    @Test
    public void testEffectiveTpmCalcs()
    {
        double reqProbLow = 0.05;
        double reqProbHigh = 0.95;

        double lowValue = TpmCalculator.calcPoissonObservedGivenProb(10, reqProbLow);
        double highValue = TpmCalculator.calcPoissonObservedGivenProb(10, reqProbHigh);

        assertEquals(4.5, lowValue, 0.1);
        assertEquals(15.0, highValue, 0.1);

        lowValue = TpmCalculator.calcPoissonObservedGivenProb(1, reqProbLow);
        highValue = TpmCalculator.calcPoissonObservedGivenProb(1, reqProbHigh);

        assertEquals(0, lowValue, 0.1);
        assertEquals(2.5, highValue, 0.1);

        lowValue = TpmCalculator.calcPoissonObservedGivenProb(2, reqProbLow);
        highValue = TpmCalculator.calcPoissonObservedGivenProb(2, reqProbHigh);

        assertEquals(0, lowValue, 0.1);
        assertEquals(4, highValue, 0.1);

        lowValue = TpmCalculator.calcPoissonObservedGivenProb(3, reqProbLow);
        highValue = TpmCalculator.calcPoissonObservedGivenProb(3, reqProbHigh);

        assertEquals(1, lowValue, 0.1);
        assertEquals(5.7, highValue, 0.1);
    }
}
