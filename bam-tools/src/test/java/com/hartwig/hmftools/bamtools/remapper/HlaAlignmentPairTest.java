package com.hartwig.hmftools.bamtools.remapper;

import static org.junit.Assert.assertEquals;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.junit.Test;

public class HlaAlignmentPairTest extends RemapperTestBase
{
    @Test
    public void interPairDistance()
    {
        BwaMemAlignment baseLeft = bwa("81,0,194358862,194358978,35,151,60,0,116,66,35S116M,116,null,5,29945439,0");
        BwaMemAlignment baseRight = bwa("161,5,29945439,29945590,0,151,60,0,151,116,151M,151,null,0,194358862,0");
        HlaAlignmentPair pair = new HlaAlignmentPair(new HlaAlignment(baseLeft), new HlaAlignment(baseRight));
        assertEquals((194358862 + 1) - (29945439 + 1), pair.interPairDistance());
    }
}
