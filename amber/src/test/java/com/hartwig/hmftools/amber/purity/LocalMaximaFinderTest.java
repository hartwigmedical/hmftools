package com.hartwig.hmftools.amber.purity;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import org.junit.Before;
import org.junit.Test;

public class LocalMaximaFinderTest
{
    private List<DummyScore> mInputs;

    @Before
    public void setup()
    {
        mInputs = new ArrayList<>();
    }

    @Test
    public void increasingTest()
    {
        addScore("0", 1.0);
        addScore("1", 2.0);
        addScore("2", 3.0);
        assertTrue(getMaxima().isEmpty());
    }

    @Test
    public void decreasingTest()
    {
        addScore("0", 10.0);
        addScore("1", 5.0);
        addScore("2", 3.0);
        assertTrue(getMaxima().isEmpty());
    }

    @Test
    public void constantTest()
    {
        addScore("0", 10.0);
        addScore("1", 10.0);
        addScore("2", 10.0);
        assertTrue(getMaxima().isEmpty());
    }

    @Test
    public void emptyTest()
    {
        assertTrue(getMaxima().isEmpty());
    }

    @Test
    public void onePointTest()
    {
        addScore("A", 10.0);
        assertTrue(getMaxima().isEmpty());
    }

    @Test
    public void plateauTest()
    {
        addScore("0", 1.0);
        addScore("1", 10.0);
        addScore("2", 10.0);
        addScore("3", 5.0);
        assertEquals(List.of(mInputs.get(1)), getMaxima());
    }

    @Test
    public void valleyTest()
    {
        addScore("0", 100.0);
        addScore("1", 10.0);
        addScore("2", 10.0);
        addScore("3", 50.0);
        assertTrue(getMaxima().isEmpty());
    }

    @Test
    public void skipInitialZeroesTest()
    {
        addScore("0", 0.000);
        addScore("1", 0.000);
        addScore("2", 0.001);
        addScore("3", 0.000);
        List<DummyScore> maxima = getMaxima();
        assertEquals(0, maxima.size());
    }

    @Test
    public void skipInitialZeroes2Test()
    {
        addScore("0", 0.000);
        addScore("1", 0.000);
        addScore("2", 0.002);
        addScore("3", 0.001);
        addScore("4", 0.000);
        List<DummyScore> maxima = getMaxima();
        assertEquals(0, maxima.size());
    }

    @Test
    public void skipInitialZeroes3Test()
    {
        addScore("0", 0.000);
        addScore("1", 0.000);
        addScore("2", 0.001);
        addScore("3", 0.002);
        addScore("4", 0.000);
        List<DummyScore> maxima = getMaxima();
        assertEquals(1, maxima.size());
        assertEquals(maxima.get(0), mInputs.get(3));
    }

    @Test
    public void typicalCaseTest()
    {
        addScore("0", 100.0);
        addScore("1", 50);
        addScore("2", 20);
        addScore("3", 18);
        addScore("4", 19);
        addScore("5", 20);
        addScore("6", 25);
        addScore("7", 30); // peak
        addScore("8", 29);
        addScore("9", 25);
        addScore("10", 21);
        addScore("11", 22.1); // false peak, will be removed
        addScore("12", 22.0); // trough
        addScore("13", 25);
        addScore("14", 30); // peak
        addScore("15", 28);
        addScore("16", 25);
        addScore("17", 30);

        List<DummyScore> maxima = getMaxima();
        assertEquals(2, maxima.size());
        assertEquals(maxima.get(0), mInputs.get(7));
        assertEquals(maxima.get(1), mInputs.get(14));
    }

    private List<DummyScore> getMaxima()
    {
        return new LocalMaximaFinder<>(mInputs).maxima();
    }

    private void addScore(String label, double score)
    {
        mInputs.add(new DummyScore(score, label));
    }

    private class DummyScore implements Score
    {
        private final double mScore;
        private final String mLabel;

        DummyScore(final double score, final String label)
        {
            mScore = score;
            mLabel = label;
        }

        @Override
        public double score()
        {
            return mScore;
        }

        @Override
        public boolean equals(final Object o)
        {
            if(o == null || getClass() != o.getClass())
            {
                return false;
            }
            final DummyScore that = (DummyScore) o;
            return Double.compare(mScore, that.mScore) == 0 && Objects.equals(mLabel, that.mLabel);
        }

        @Override
        public int hashCode()
        {
            return Objects.hash(mScore, mLabel);
        }

        @Override
        public String toString()
        {
            return mLabel;
        }
    }
}

