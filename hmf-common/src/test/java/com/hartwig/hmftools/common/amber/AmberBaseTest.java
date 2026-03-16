package com.hartwig.hmftools.common.amber;

import static com.hartwig.hmftools.common.amber.AmberBase.A;
import static com.hartwig.hmftools.common.amber.AmberBase.C;
import static com.hartwig.hmftools.common.amber.AmberBase.G;
import static com.hartwig.hmftools.common.amber.AmberBase.N;
import static com.hartwig.hmftools.common.amber.AmberBase.T;

import org.junit.Assert;
import org.junit.Test;

public class AmberBaseTest
{
    @Test
    public void complementTest()
    {
        Assert.assertEquals(T, A.complement());
        Assert.assertEquals(A, T.complement());
        Assert.assertEquals(G, C.complement());
        Assert.assertEquals(C, G.complement());
        Assert.assertEquals(N, N.complement());
    }
}
