package com.hartwig.hmftools.cobalt.count;

import static org.mockito.ArgumentMatchers.any;
import static org.mockito.ArgumentMatchers.anyInt;
import static org.mockito.ArgumentMatchers.anyString;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.never;
import static org.mockito.Mockito.reset;
import static org.mockito.Mockito.verify;

import static htsjdk.samtools.util.SequenceUtil.A;
import static htsjdk.samtools.util.SequenceUtil.C;
import static htsjdk.samtools.util.SequenceUtil.G;
import static htsjdk.samtools.util.SequenceUtil.T;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class BamReadCounterTest
{
    @Test
    public void testAccumulateAlignmentBlock()
    {
        final String CHROMOSOME = "chr1";

        final ReadDepthAccumulator readDepthAccumulator = mock(ReadDepthAccumulator.class);
        final byte[] readBases = new byte[210];

        // create some base patterns
        for(int i = 0; i < readBases.length; ++i)
        {
            switch(i % 7)
            {
                case 0:
                case 3:
                    readBases[i] = A; break;
                case 1: readBases[i] = G; break;
                case 2:
                case 5:
                    readBases[i] = C; break;
                default: readBases[i] = T; break;
            }
        }
        
        final ChrBaseRegion region = new ChrBaseRegion(CHROMOSOME, 1001, 2000);

        BamReadCounter.accumulateAlignmentBlock(region, readDepthAccumulator,
                11, 1101, 200, readBases);
        // check that it invoked the correct parameters
        verify(readDepthAccumulator).addReadAlignmentToCounts(CHROMOSOME, 1101, 200, readBases, 10);

        // now test that first 100 bases are before this region
        BamReadCounter.accumulateAlignmentBlock(region, readDepthAccumulator,
                11, 901, 200, readBases);
        // check that it invoked the correct parameters
        verify(readDepthAccumulator).addReadAlignmentToCounts(CHROMOSOME, 1001, 100, readBases, 110);

        // now test that last 100 bases are after this region
        BamReadCounter.accumulateAlignmentBlock(region, readDepthAccumulator,
                11, 1901, 200, readBases);
        // check that it invoked the correct parameters
        verify(readDepthAccumulator).addReadAlignmentToCounts(CHROMOSOME, 1901, 100, readBases, 10);

        reset(readDepthAccumulator);

        // test a block that is before the region, check that no call to addReadAlignmentToCounts has been invoked
        BamReadCounter.accumulateAlignmentBlock(region, readDepthAccumulator,
                11, 500, 200, readBases);
        verify(readDepthAccumulator, never()).addReadAlignmentToCounts(anyString(), anyInt(), anyInt(), any(byte[].class), anyInt());

        // test a block that is after the region, check that no call to addReadAlignmentToCounts has been invoked
        BamReadCounter.accumulateAlignmentBlock(region, readDepthAccumulator,
                11, 2500, 200, readBases);
        verify(readDepthAccumulator, never()).addReadAlignmentToCounts(anyString(), anyInt(), anyInt(), any(byte[].class), anyInt());
    }
}
