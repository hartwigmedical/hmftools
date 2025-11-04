package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertSame;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import com.hartwig.hmftools.cobalt.consolidation.NoOpConsolidator;
import com.hartwig.hmftools.cobalt.consolidation.ResultsConsolidator;
import com.hartwig.hmftools.cobalt.normalisers.DoNothingNormaliser;
import com.hartwig.hmftools.cobalt.normalisers.ReadDepthStatisticsNormaliser;
import com.hartwig.hmftools.cobalt.normalisers.ResultsNormaliser;
import com.hartwig.hmftools.cobalt.targeted.CobaltScope;

import org.junit.Assert;
import org.junit.Test;

public class TumorCalculationTest
{
    private TumorCalculation Calculation;

    @Test
    public void createFinalNormaliser()
    {
        ResultsNormaliser theLast = mock(ResultsNormaliser.class);
        CobaltScope scope = mock(CobaltScope.class);
        when(scope.finalNormaliser()).thenReturn(theLast);
        Calculation = new TumorCalculation(mock(WindowStatuses.class), scope, V38);
        assertEquals(theLast, Calculation.createFinalNormaliser());
    }

    @Test
    public void consolidator()
    {
        ResultsConsolidator consolidator = new NoOpConsolidator();
        ReadDepthStatisticsNormaliser readDepthStatisticsNormaliser = mock(ReadDepthStatisticsNormaliser.class);
        when(readDepthStatisticsNormaliser.readDepthMedian()).thenReturn(4.0);
        CobaltScope scope = mock(CobaltScope.class);
        when(scope.medianByMeanNormaliser()).thenReturn(readDepthStatisticsNormaliser);
        when(scope.resultsConsolidator(4.0)).thenReturn(consolidator);
        Calculation = new TumorCalculation(mock(WindowStatuses.class), scope, V38);
        assertEquals(consolidator, Calculation.consolidator());
        // Subsequent calls should return the same instance.
        assertSame(consolidator, Calculation.consolidator());
    }

    @Test
    public void createReadDepthsNormaliser()
    {
        ReadDepthStatisticsNormaliser theMeanOne = mock(ReadDepthStatisticsNormaliser.class);
        CobaltScope scope = mock(CobaltScope.class);
        when(scope.medianByMeanNormaliser()).thenReturn(theMeanOne);
        Calculation = new TumorCalculation(mock(WindowStatuses.class), scope, V38);
        assertEquals(theMeanOne, Calculation.createReadDepthsNormaliser());
    }

    @Test
    public void createMegaBaseScaleNormaliser()
    {
        Calculation = new TumorCalculation(mock(WindowStatuses.class), mock(CobaltScope.class), V38);
        Assert.assertTrue(Calculation.createMegaBaseScaleNormaliser(V38) instanceof DoNothingNormaliser);
    }
}
