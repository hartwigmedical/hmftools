package com.hartwig.hmftools.cobalt.calculations;

import static org.mockito.Mockito.when;

import com.hartwig.hmftools.cobalt.targeted.CobaltScope;

import org.junit.Assert;
import org.junit.Test;
import org.mockito.Mockito;

public class TumorBamCalculationTest
{
    private TumorBamCalculation Calculation;

    @Test
    public void createFinalNormaliser()
    {
        ResultsNormaliser theLast = Mockito.mock(ResultsNormaliser.class);
        CobaltScope scope = Mockito.mock(CobaltScope.class);
        when(scope.finalNormaliser()).thenReturn(theLast);
        Calculation = new TumorBamCalculation(Mockito.mock(GenomeFilter.class), scope);
        Assert.assertEquals(theLast, Calculation.createReadDepthsNormaliser());
    }

    @Test
    public void createMegaBaseScaleNormaliser()
    {
        Calculation = new TumorBamCalculation(Mockito.mock(GenomeFilter.class), Mockito.mock(CobaltScope.class));
        Assert.assertTrue(Calculation.createMegaBaseScaleNormaliser() instanceof DoNothingNormaliser);
    }
}
