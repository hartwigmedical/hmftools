package com.hartwig.hmftools.neo;

import static com.hartwig.hmftools.neo.cohort.AlleleCoverage.EXPECTED_ALLELE_COUNT;

import static junit.framework.TestCase.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.neo.cohort.PeptidePredictionData;

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
}
