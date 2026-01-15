package com.hartwig.hmftools.common.purple;

import static com.hartwig.hmftools.common.purple.ChromosomeArm.CENTROMERE;
import static com.hartwig.hmftools.common.purple.ChromosomeArm.P_ARM;
import static com.hartwig.hmftools.common.purple.ChromosomeArm.Q_ARM;
import static com.hartwig.hmftools.common.purple.ChromosomeArm.UNKNOWN;

import org.junit.Assert;
import org.junit.Test;

public class ChromosomeArmTest
{
    @Test
    public void asStrTest()
    {
        Assert.assertEquals("P", ChromosomeArm.asStr(P_ARM));
        Assert.assertEquals("Q", ChromosomeArm.asStr(Q_ARM));
        Assert.assertEquals("CENTROMERE", ChromosomeArm.asStr(CENTROMERE));
        Assert.assertEquals("UNKNOWN", ChromosomeArm.asStr(UNKNOWN));
    }

    @Test
    public void fromStrTest()
    {
        Assert.assertEquals(P_ARM, ChromosomeArm.fromString("P"));
        Assert.assertEquals(ChromosomeArm.Q_ARM, ChromosomeArm.fromString("Q"));
        Assert.assertEquals(ChromosomeArm.CENTROMERE, ChromosomeArm.fromString("CENTROMERE"));
        Assert.assertEquals(UNKNOWN, ChromosomeArm.fromString("UNKNOWN"));
    }
}
