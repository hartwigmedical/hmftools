package com.hartwig.hmftools.amber.contamination;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.apache.commons.lang3.tuple.Pair;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class SearchGridTest
{

    @Test
    public void searchValuesAndStepsTest()
    {
        final double delta = 0.00001;
        final SearchGrid searchGrid = new SearchGrid();
        final List<Pair<Double, Double>> searchValues = searchGrid.searchValuesAndSteps();
        assertEquals(0.0050, searchValues.get(0).getLeft(), delta);
        assertEquals(0.00105, searchValues.get(0).getRight(), delta);
        assertEquals(0.0060, searchValues.get(1).getLeft(), delta);
        assertEquals(0.0011, searchValues.get(1).getRight(), delta);
        assertEquals(0.007, searchValues.get(2).getLeft(), delta);
        assertEquals(0.00116, searchValues.get(2).getRight(), delta);
        assertEquals(0.341, searchValues.get(searchValues.size() - 1).getLeft(), delta);
        assertEquals(0.324, searchValues.get(searchValues.size() - 2).getLeft(), delta);
        assertEquals(1.05,
                searchValues.get(searchValues.size() - 1).getRight() / searchValues.get(searchValues.size() - 2).getRight(), delta);
        assertTrue(searchValues.size() >= 60);
    }

    @Test
    public void searchValuesTest()
    {
        final SearchGrid searchGrid = new SearchGrid();
        final List<Double> searchValues = searchGrid.searchValues();
        assertEquals(0.005, searchValues.get(0), 0.0001);
        assertEquals(0.006, searchValues.get(1), 0.0001);
        assertEquals(0.007, searchValues.get(2), 0.0001);
        assertEquals(0.341, searchValues.get(searchValues.size() - 1), 0.0001);
        assertEquals(0.324, searchValues.get(searchValues.size() - 2), 0.0001);
        assertTrue(searchValues.size() >= 60);
        checkIncreasingByAtLeast(searchValues);
    }

    @Test
    public void valueScoreComparisonTest()
    {
        SearchGrid.ValueScore vs1 = new SearchGrid.ValueScore(0.05, 100);
        SearchGrid.ValueScore vs2 = new SearchGrid.ValueScore(0.05, 100);
        assertEquals(0, vs1.compareTo(vs2));
        SearchGrid.ValueScore vs3 = new SearchGrid.ValueScore(0.05, 101);
        assertTrue(vs1.compareTo(vs3) < 0);
        SearchGrid.ValueScore vs4 = new SearchGrid.ValueScore(0.04, 100);
        assertTrue(vs1.compareTo(vs4) > 0);
        SearchGrid.ValueScore vs5 = new SearchGrid.ValueScore(0.04, 101);
        assertTrue(vs1.compareTo(vs5) < 0);
        SearchGrid.ValueScore vs6 = new SearchGrid.ValueScore(0.04, 99);
        assertTrue(vs1.compareTo(vs6) > 0);
        SearchGrid.ValueScore vs7 = new SearchGrid.ValueScore(0.06, 99);
        assertTrue(vs1.compareTo(vs7) > 0);
    }

    private void checkIncreasingByAtLeast(List<Double> values)
    {
        double previous = -1.0;
        for(Double v : values)
        {
            assertTrue(v - previous >= 0.001);
        }
    }
}
