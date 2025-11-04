package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import com.hartwig.hmftools.cobalt.consolidation.NoOpConsolidator;
import com.hartwig.hmftools.cobalt.normalisers.DiploidNormaliser;
import com.hartwig.hmftools.cobalt.normalisers.ReadDepthStatisticsNormaliser;
import com.hartwig.hmftools.cobalt.targeted.CobaltScope;

import org.junit.Assert;
import org.junit.Test;

public class ReferenceCalculationTest
{
    private ReferenceCalculation Calculation;

    @Test
    public void createReadDepthsNormaliser()
    {
        ReadDepthStatisticsNormaliser theMeanOne = mock(ReadDepthStatisticsNormaliser.class);
        CobaltScope scope = mock(CobaltScope.class);
        when(scope.medianByMeanNormaliser()).thenReturn(theMeanOne);
        Calculation = new ReferenceCalculation(mock(WindowStatuses.class), scope, V38, new NoOpConsolidator());
        Assert.assertEquals(theMeanOne, Calculation.createReadDepthsNormaliser());
    }

    @Test
    public void createMegaBaseScaleNormaliser()
    {
        Calculation = new ReferenceCalculation(mock(WindowStatuses.class), mock(CobaltScope.class), V38, new NoOpConsolidator());
        Assert.assertTrue(Calculation.createMegaBaseScaleNormaliser(V38) instanceof DiploidNormaliser);
    }
}
