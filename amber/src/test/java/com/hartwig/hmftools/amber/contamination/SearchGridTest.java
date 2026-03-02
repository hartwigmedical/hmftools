package com.hartwig.hmftools.amber.contamination;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class SearchGridTest
{
    @Test
    public void closeToZeroTest()
    {
        final SearchGrid searchGrid = new SearchGrid();
        final List<Double> searchValues = searchGrid.searchValues();
        assertEquals(0.01, searchValues.get(0), 0.0001);
        assertEquals(0.06, searchValues.get(searchValues.size() - 1), 0.0001);
        assertTrue(searchValues.size() > 25);
        checkIncreasingByAtLeast(searchValues);
    }

    @Test
    public void searchValuesTest()
    {
        final SearchGrid searchGrid = new SearchGrid();
        final List<Double> searchValues = searchGrid.searchValues();
        assertEquals(0.05, searchValues.get(0), 0.0001);
        assertEquals(0.053, searchValues.get(1), 0.0001);
        assertEquals(0.055, searchValues.get(2), 0.0001);
        assertEquals(0.335, searchValues.get(searchValues.size() - 1), 0.0001);
        assertEquals(0.319, searchValues.get(searchValues.size() - 2), 0.0001);
        assertTrue(searchValues.size() > 50);
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
