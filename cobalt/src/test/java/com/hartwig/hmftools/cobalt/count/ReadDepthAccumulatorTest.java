package com.hartwig.hmftools.cobalt.count;

import static com.hartwig.hmftools.cobalt.CobaltTestUtils.EPSILON;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import static htsjdk.samtools.util.SequenceUtil.A;
import static htsjdk.samtools.util.SequenceUtil.C;
import static htsjdk.samtools.util.SequenceUtil.G;
import static htsjdk.samtools.util.SequenceUtil.T;

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
        List<DepthReading> readDepths = readDepthCounter.getChromosomeReadDepths(CHROMOSOME);

        assertNotNull(readDepths);
        assertEquals(2, readDepths.size());

        byte[] bases = new byte[1000];

        // we fill first 400 bases with GC, and next 400 base with AT, then last 200 bases with GC
        for(int i = 0; i < bases.length; ++i)
        {
            if(i < 400 || i >= 800)
            {
                if (i % 2 == 0)
                    bases[i] = G;
                else
                    bases[i] = C;
            }
            else if (i % 2 == 0)
                bases[i] = A;
            else
                bases[i] = T;
        }

        // add some read data
        readDepthCounter.addReadAlignmentToCounts(CHROMOSOME, 501, 1000, bases, 0);

        readDepths = readDepthCounter.getChromosomeReadDepths(CHROMOSOME);
        assertNotNull(readDepths);
        assertEquals(2, readDepths.size());
        DepthReading readDepth = readDepths.get(0);
        assertEquals(1, readDepth.StartPosition);
        assertEquals(0.5, readDepth.ReadDepth, EPSILON);
        // gc percent should be 0.8 as first 400 bases were GC, and next 100 bases were AT
        assertEquals(0.8, readDepth.ReadGcContent, EPSILON);
        readDepth = readDepths.get(1);
        assertEquals(1001, readDepth.StartPosition);
        assertEquals(0.5, readDepth.ReadDepth, EPSILON);
        // gc percent should be 0.4 as first 300 bases were AT, and next 200 bases were GC
        assertEquals(0.4, readDepth.ReadGcContent, EPSILON);

        // add one more read that only covers the first window

        // fill it with GAGAGAGAGAGA...
        for(int i = 0; i < bases.length; ++i)
        {
            if(i % 2 == 0)
                bases[i] = G;
            else
                bases[i] = A;
        }

        readDepthCounter.addReadAlignmentToCounts(CHROMOSOME, 1, 1000, bases, 0);
        readDepths = readDepthCounter.getChromosomeReadDepths(CHROMOSOME);
        assertNotNull(readDepths);
        assertEquals(2, readDepths.size());
        readDepth = readDepths.get(0);
        assertEquals(1, readDepth.StartPosition);
        assertEquals(1.5, readDepth.ReadDepth, EPSILON);

        // first read has 0.8 gc with 500 bases covered, second read has 0.5 gc with 1000 bases
        // together it becomes 0.6
        assertEquals(0.6, readDepth.ReadGcContent, EPSILON);

        // second read should have no change
        readDepth = readDepths.get(1);
        assertEquals(1001, readDepth.StartPosition);
        assertEquals(0.5, readDepth.ReadDepth, EPSILON);
        assertEquals(0.4, readDepth.ReadGcContent, EPSILON);
    }
}
