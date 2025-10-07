package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static org.mockito.Mockito.when;

import com.hartwig.hmftools.cobalt.targeted.CobaltScope;

import org.junit.Assert;
import org.junit.Test;
import org.mockito.Mockito;

public class ReferenceBamCalculationTest
{
    private ReferenceBamCalculation Calculation;

    @Test
    public void createReadDepthsNormaliser()
    {
        ReadDepthStatisticsNormaliser theMeanOne = Mockito.mock(ReadDepthStatisticsNormaliser.class);
        CobaltScope scope = Mockito.mock(CobaltScope.class);
        when(scope.medianByMeanNormaliser()).thenReturn(theMeanOne);
        Calculation = new ReferenceBamCalculation(Mockito.mock(GenomeFilter.class), scope, V38);
        Assert.assertEquals(theMeanOne, Calculation.createReadDepthsNormaliser());
    }

    @Test
    public void createMegaBaseScaleNormaliser()
    {
        Calculation = new ReferenceBamCalculation(Mockito.mock(GenomeFilter.class), Mockito.mock(CobaltScope.class), V38);
        Assert.assertTrue(Calculation.createMegaBaseScaleNormaliser(Mockito.any()) instanceof DiploidNormaliser);
    }
}
