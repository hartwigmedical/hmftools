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
    private List<DummyScore> Inputs;

    @Before
    public void setup()
    {
        Inputs = new ArrayList<>();
    }

    @Test
    public void increasingTest()
    {
        addScore("A", 1.0);
        addScore("B", 2.0);
        addScore("C", 3.0);
        assertTrue(getMaxima().isEmpty());
    }

    @Test
    public void decreasingTest()
    {
        addScore("A", 10.0);
        addScore("B", 5.0);
        addScore("C", 3.0);
        assertTrue(getMaxima().isEmpty());
    }

    @Test
    public void constantTest()
    {
        addScore("A", 10.0);
        addScore("B", 10.0);
        addScore("C", 10.0);
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
        addScore("A", 1.0);
        addScore("B", 10.0);
        addScore("C", 10.0);
        addScore("D", 5.0);
        assertEquals(List.of(Inputs.get(2)), getMaxima());
    }

    @Test
    public void valleyTest()
    {
        addScore("A", 100.0);
        addScore("B", 10.0);
        addScore("C", 10.0);
        addScore("D", 50.0);
        assertTrue(getMaxima().isEmpty());
    }

    @Test
    public void skipInitialZeroesTest()
    {
        addScore("A", 0.000);
        addScore("B", 0.000);
        addScore("C", 0.001);
        addScore("D", 0.000);
        List<DummyScore> maxima = getMaxima();
        assertEquals(0, maxima.size());
    }

    @Test
    public void skipInitialZeroes2Test()
    {
        addScore("A", 0.000);
        addScore("B", 0.000);
        addScore("C", 0.002);
        addScore("D", 0.001);
        addScore("E", 0.000);
        List<DummyScore> maxima = getMaxima();
        assertEquals(0, maxima.size());
    }

    @Test
    public void skipInitialZeroes3Test()
    {
        addScore("A", 0.000);
        addScore("B", 0.000);
        addScore("C", 0.001);
        addScore("D", 0.002);
        addScore("E", 0.000);
        List<DummyScore> maxima = getMaxima();
        assertEquals(1, maxima.size());
        assertEquals(maxima.get(0), Inputs.get(3));
    }

    @Test
    public void typicalCaseTest()
    {
        addScore("A", 100.0);
        addScore("B", 50);
        addScore("C", 20);
        addScore("D", 18);
        addScore("E", 19);
        addScore("F", 20);
        addScore("G", 25);
        addScore("H", 30); // *
        addScore("I", 29);
        addScore("J", 25);
        addScore("K", 22);
        addScore("L", 18);
        addScore("M", 20);
        addScore("N", 25);
        addScore("O", 30); // *
        addScore("P", 28);
        addScore("Q", 25);
        addScore("R", 30); // ? TODO
        List<DummyScore> maxima = getMaxima();
        assertEquals(2, maxima.size());
        assertEquals(maxima.get(0), Inputs.get(7));
        assertEquals(maxima.get(1), Inputs.get(14));
    }

    private List<DummyScore> getMaxima()
    {
        return new LocalMaximaFinder<>(Inputs).maxima();
    }

    private void addScore(String label, double score)
    {
        Inputs.add(new DummyScore(score, label));
    }
}

class DummyScore implements Score
{
    private final double mScore;
    private final String Label;

    DummyScore(final double mScore, final String label)
    {
        this.mScore = mScore;
        Label = label;
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
        return Double.compare(mScore, that.mScore) == 0 && Objects.equals(Label, that.Label);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mScore, Label);
    }

    @Override
    public String toString()
    {
        return Label;
    }
}
