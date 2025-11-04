package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.cobalt.consolidation.NoOpConsolidator;
import com.hartwig.hmftools.cobalt.consolidation.ResultsConsolidator;
import com.hartwig.hmftools.cobalt.normalisers.ReadDepthStatisticsNormaliser;
import com.hartwig.hmftools.cobalt.normalisers.ResultsNormaliser;
import com.hartwig.hmftools.cobalt.targeted.CobaltScope;
import com.hartwig.hmftools.common.cobalt.GcMedianReadDepth;
import com.hartwig.hmftools.common.genome.gc.ImmutableGCBucket;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.junit.Before;
import org.junit.Test;
import org.mockito.Mockito;

public class BamCalculationTest extends CalculationsTestBase
{
    RecordingNormaliser DepthNormaliser;
    RecordingNormaliser MbScaleNormaliser;
    RecordingNormaliser FinalNormaliser;
    BC calculation;

    @Before
    public void setup()
    {
        DepthNormaliser = new RecordingNormaliser(1.0 / 2.0);
        MbScaleNormaliser = new RecordingNormaliser(1.0 / 3.0);
        FinalNormaliser = new RecordingNormaliser(1.0 / 5.0);
        final WindowStatuses genomeFilter = Mockito.mock(WindowStatuses.class);
        Mockito.when(genomeFilter.exclude(Mockito.any(), Mockito.any())).thenReturn(false);
        final CobaltScope scope = Mockito.mock(CobaltScope.class);
        Mockito.when(scope.onTarget(Mockito.any(), Mockito.any(Integer.class))).thenReturn(true);
        Mockito.when(scope.enrichmentQuotient(Mockito.any(), Mockito.anyInt())).thenReturn(1.0);
        calculation = new BC(genomeFilter, scope, V38);
    }

    @Test
    public void normalisationOrderTest()
    {
        calculation.addReading(_1, dr(_1, 1, 1, 0.40));

        calculation.calculateRatios();

        assertEquals(1, DepthNormaliser.RecordedValues.size());
        // Input is 1.0. The GC depth is 1/3 because of GC bucket smoothing (there's only one bucket).
        // So the value passed to the depth normaliser should be 1.0/1/3 = 3.0.
        assertEquals(3.0, DepthNormaliser.RecordedValues.get(0), 0.0001);

        // The input to the mb scale normaliser should be the output of the depth normaliser.
        assertEquals(1, MbScaleNormaliser.RecordedValues.size());
        assertEquals(3.0 * 2.0, MbScaleNormaliser.RecordedValues.get(0), 0.0001);

        // The input to the final normaliser should be the output of the mb scale normaliser.
        assertEquals(1, FinalNormaliser.RecordedValues.size());
        assertEquals(3.0 * 3.0 * 2.0, FinalNormaliser.RecordedValues.get(0), 0.0001);
    }

    @Test
    public void medianReadDepthsTest()
    {
        calculation.addReading(_1, dr(_1, 1, 1, 0.40));
        calculation.calculateRatios();
        GcMedianReadDepth gcMedianReadDepth = calculation.medianReadDepths();
        assertEquals(1.23, gcMedianReadDepth.meanReadDepth(), 0.001);
        assertEquals(2.45, gcMedianReadDepth.medianReadDepth(), 0.001);
        assertEquals(1.0 / 3.0, gcMedianReadDepth.medianReadDepth(new ImmutableGCBucket(40)), 0.001);
    }

    class BC extends BamCalculation
    {
        public BC(final WindowStatuses mGenomeFilter, final CobaltScope scope, final RefGenomeVersion version)
        {
            super(mGenomeFilter, scope, version);
        }

        @Override
        ReadDepthStatisticsNormaliser createReadDepthsNormaliser()
        {
            return DepthNormaliser;
        }

        @Override
        ResultsNormaliser createMegaBaseScaleNormaliser(final RefGenomeVersion version)
        {
            return MbScaleNormaliser;
        }

        @Override
        ResultsNormaliser createFinalNormaliser()
        {
            return FinalNormaliser;
        }

        @Override
        ResultsConsolidator consolidator()
        {
            return new NoOpConsolidator();
        }
    }

    static class RecordingNormaliser extends ReadDepthStatisticsNormaliser
    {
        final double Factor;
        final List<Double> NormalisedValues = new ArrayList<>();
        final List<Double> RecordedValues = new ArrayList<>();

        RecordingNormaliser(final double factor)
        {
            this.Factor = factor;
        }

        @Override
        public void recordValue(final BamRatio bamRatio)
        {
            RecordedValues.add(bamRatio.ratio());
        }

        @Override
        public void normalise(final BamRatio bamRatio)
        {
            NormalisedValues.add(bamRatio.ratio());
            bamRatio.normaliseByMean(Factor);
        }

        @Override
        public double readDepthMean()
        {
            return 1.23;
        }

        @Override
        public double readDepthMedian()
        {
            return 2.45;
        }
    }
}
