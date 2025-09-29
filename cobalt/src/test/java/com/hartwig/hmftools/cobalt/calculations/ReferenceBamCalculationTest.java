package com.hartwig.hmftools.cobalt.calculations;

import com.hartwig.hmftools.cobalt.targeted.CobaltScope;

import org.junit.Assert;
import org.junit.Test;
import org.mockito.Mockito;

public class ReferenceBamCalculationTest
{
    private ReferenceBamCalculation Calculation;

    @Test
    public void finalNormaliser()
    {
        Calculation = new ReferenceBamCalculation(Mockito.mock(GenomeFilter.class), Mockito.mock(CobaltScope.class));
        Assert.assertTrue(Calculation.finalMeanNormaliser() instanceof DoNothingNormaliser);
    }

    @Test
    public void diploidNormaliser()
    {
        Calculation = new ReferenceBamCalculation(Mockito.mock(GenomeFilter.class), Mockito.mock(CobaltScope.class));
        Assert.assertTrue(Calculation.diploidNormaliser() instanceof DiploidNormaliser);
    }
}
