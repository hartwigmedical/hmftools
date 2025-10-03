package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import com.hartwig.hmftools.cobalt.targeted.CobaltScope;

import org.junit.Assert;
import org.junit.Test;
import org.mockito.Mockito;

public class ReferenceBamCalculationTest
{
    private ReferenceBamCalculation Calculation;

    @Test
    public void createFinalNormaliser()
    {
        Calculation = new ReferenceBamCalculation(Mockito.mock(GenomeFilter.class), Mockito.mock(CobaltScope.class), V38);
        Assert.assertTrue(Calculation.createReadDepthsNormaliser() instanceof ReadDepthStatisticsNormaliser);
    }

    @Test
    public void createMegaBaseScaleNormaliser()
    {
        Calculation = new ReferenceBamCalculation(Mockito.mock(GenomeFilter.class), Mockito.mock(CobaltScope.class), V38);
        Assert.assertTrue(Calculation.createMegaBaseScaleNormaliser(Mockito.any()) instanceof DiploidNormaliser);
    }
}
