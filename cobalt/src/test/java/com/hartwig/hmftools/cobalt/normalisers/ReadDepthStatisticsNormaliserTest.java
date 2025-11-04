package com.hartwig.hmftools.cobalt.normalisers;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._X;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._Y;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.cobalt.calculations.BamRatio;
import com.hartwig.hmftools.cobalt.calculations.CalculationsTestBase;

import org.junit.Before;
import org.junit.Test;

public class ReadDepthStatisticsNormaliserTest extends CalculationsTestBase
{
    private ReadDepthStatisticsNormaliser mNormaliser;

    @Before
    public void setup()
    {
        mNormaliser = new ReadDepthStatisticsNormaliser();
    }

    @Test
    public void normaliseTest()
    {
        mNormaliser.recordValue(br(_1, 1, 10, 0.5, true));
        mNormaliser.recordValue(br(_1, 101, 10, 0.5, true));
        mNormaliser.recordValue(br(_1, 201, 100, 0.5, true));

        mNormaliser.dataCollectionFinished();

        // Mean = 40, median = 10. Expected factor = 0.25
        BamRatio ratio =  br(_1, 1, 1.0, 0.5, true);
        assertEquals(1.0, ratio.ratio(), 0.0001); // sanity check
        mNormaliser.normalise(ratio);
        assertEquals(0.25, ratio.ratio(), 0.0001);
    }

    @Test
    public void valuesWithExtremeGcAreIgnored()
    {
        mNormaliser.recordValue(br(_1, 1, 100, 0.23, true));
        mNormaliser.recordValue(br(_1, 6, 100, 0.23999, true));
        mNormaliser.recordValue(br(_1, 11, 10, 0.24, true));
        mNormaliser.recordValue(br(_1, 21, 10, 0.25, true));
        mNormaliser.recordValue(br(_1, 31, 10, 0.26, true));
        mNormaliser.recordValue(br(_1, 41, 10, 0.34, true));
        mNormaliser.recordValue(br(_1, 51, 10, 0.54, true));
        mNormaliser.recordValue(br(_1, 61, 10, 0.67, true));
        mNormaliser.recordValue(br(_1, 71, 10, 0.68, true));
        mNormaliser.recordValue(br(_1, 71, 1000, 0.680001, true));
        mNormaliser.recordValue(br(_1, 81, 1000, 0.69, true));

        mNormaliser.dataCollectionFinished();

        // In-range depths are all 10, so norm does nothing.
        BamRatio ratio =  br(_1, 1, 1.0, 0.5, true);
        mNormaliser.normalise(ratio);
        assertEquals(1.0, ratio.ratio(), 0.0001);
    }

    @Test
    public void sexChromosomeValuesAreIgnored()
    {
        mNormaliser.recordValue(br(_1, 11, 10, 0.44, true));
        mNormaliser.recordValue(br(_1, 21, 10, 0.45, true));
        mNormaliser.recordValue(br(_1, 31, 10, 0.46, true));
        mNormaliser.recordValue(br(_1, 41, 10, 0.34, true));
        mNormaliser.recordValue(br(_1, 51, 10, 0.54, true));
        mNormaliser.recordValue(br(_X, 61, 1000, 0.67, true));
        mNormaliser.recordValue(br(_Y, 71, 10000, 0.48, true));

        mNormaliser.dataCollectionFinished();

        // autosome depths are all 10, so norm does nothing.
        BamRatio ratio =  br(_1, 1, 1.0, 0.5, true);
        mNormaliser.normalise(ratio);
        assertEquals(1.0, ratio.ratio(), 0.0001);
    }

    @Test
    public void nonPositiveValuesIgnored()
    {
        mNormaliser.recordValue(br(_1, 11, 10, 0.44, true));
        mNormaliser.recordValue(br(_1, 21, 10, 0.45, true));
        mNormaliser.recordValue(br(_1, 31, 10, 0.46, true));
        mNormaliser.recordValue(br(_1, 41, 0.0, 0.34, true));
        mNormaliser.recordValue(br(_1, 51, -1.0, 0.54, true));
        mNormaliser.recordValue(br(_1, 61, -0.0001, 0.54, true));

        mNormaliser.dataCollectionFinished();

        // positive depths are all 10, so norm does nothing.
        BamRatio ratio =  br(_1, 1, 1.0, 0.5, true);
        mNormaliser.normalise(ratio);
        assertEquals(1.0, ratio.ratio(), 0.0001);
    }
}
