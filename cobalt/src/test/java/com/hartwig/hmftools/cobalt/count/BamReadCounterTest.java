package com.hartwig.hmftools.cobalt.count;

import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.verify;

import static htsjdk.samtools.util.SequenceUtil.A;
import static htsjdk.samtools.util.SequenceUtil.C;
import static htsjdk.samtools.util.SequenceUtil.G;
import static htsjdk.samtools.util.SequenceUtil.T;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class BamReadCounterTest
{
    @Test
    public void testAccumulateAlignmentBlock()
    {
        String CHROMOSOME = "chr1";

        ReadDepthAccumulator readDepthCounter = mock(ReadDepthAccumulator.class);
        byte[] bases = new byte[210];

        for(int i = 0; i < bases.length; ++i)
        {
            switch(i % 7)
            {
                case 0:
                case 3:
                    bases[i] = A; break;
                case 1: bases[i] = G; break;
                case 2:
                case 5:
                    bases[i] = C; break;
                default: bases[i] = T; break;
            }
        }

        BamReadCounter.accumulateAlignmentBlock(new ChrBaseRegion(CHROMOSOME, 1001, 2000), readDepthCounter,
                11, 1101, 200, bases);

        // check that it invoked the correct parameters
        verify(readDepthCounter).addReadAlignmentToCounts(CHROMOSOME, 1101, 200, bases, 10);

        // now test that first 100 bases are before this region
        BamReadCounter.accumulateAlignmentBlock(new ChrBaseRegion(CHROMOSOME, 1001, 2000), readDepthCounter,
                11, 901, 200, bases);

        // check that it invoked the correct parameters
        verify(readDepthCounter).addReadAlignmentToCounts(CHROMOSOME, 1001, 100, bases, 110);

        // now test that last 100 bases are after this region
        BamReadCounter.accumulateAlignmentBlock(new ChrBaseRegion(CHROMOSOME, 1001, 2000), readDepthCounter,
                11, 1901, 200, bases);

        // check that it invoked the correct parameters
        verify(readDepthCounter).addReadAlignmentToCounts(CHROMOSOME, 1901, 100, bases, 10);
    }
}
