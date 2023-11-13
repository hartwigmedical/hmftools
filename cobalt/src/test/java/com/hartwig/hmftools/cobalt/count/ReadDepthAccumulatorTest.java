package com.hartwig.hmftools.cobalt.count;

import static com.hartwig.hmftools.cobalt.CobaltTestUtils.EPSILON;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.util.List;

import org.junit.Test;

public class ReadDepthAccumulatorTest
{
    @Test
    public void testWindowReadDepth()
    {
        String CHROMOSOME = "chr1";

        ReadDepthAccumulator readDepthCounter = new ReadDepthAccumulator(1000);
        readDepthCounter.addChromosome(CHROMOSOME, 2000);

        // test that it gets the correct windows
        List<ReadDepth> readDepths = readDepthCounter.getChromosomeReadDepths(CHROMOSOME);

        assertNotNull(readDepths);
        assertEquals(2, readDepths.size());

        // add some read data
        readDepthCounter.addReadAlignmentToCounts(CHROMOSOME, 501, 1500);

        readDepths = readDepthCounter.getChromosomeReadDepths(CHROMOSOME);
        assertNotNull(readDepths);
        assertEquals(2, readDepths.size());
        ReadDepth readDepth = readDepths.get(0);
        assertEquals(1, readDepth.startPosition);
        assertEquals(0.5, readDepth.readDepth, EPSILON);
        readDepth = readDepths.get(1);
        assertEquals(1001, readDepth.startPosition);
        assertEquals(0.5, readDepth.readDepth, EPSILON);

        // add one more read that only covers the first window
        readDepthCounter.addReadAlignmentToCounts(CHROMOSOME, 1, 1000);
        readDepths = readDepthCounter.getChromosomeReadDepths(CHROMOSOME);
        assertNotNull(readDepths);
        assertEquals(2, readDepths.size());
        readDepth = readDepths.get(0);
        assertEquals(1, readDepth.startPosition);
        assertEquals(1.5, readDepth.readDepth, EPSILON);
        readDepth = readDepths.get(1);
        assertEquals(1001, readDepth.startPosition);
        assertEquals(0.5, readDepth.readDepth, EPSILON);
    }
}
