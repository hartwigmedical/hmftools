package com.hartwig.hmftools.amber.purity;

import static com.hartwig.hmftools.amber.purity.CanonicalSnvType.C_A;
import static com.hartwig.hmftools.amber.purity.CanonicalSnvType.C_G;
import static com.hartwig.hmftools.amber.purity.CanonicalSnvType.C_T;
import static com.hartwig.hmftools.amber.purity.CanonicalSnvType.T_A;
import static com.hartwig.hmftools.amber.purity.CanonicalSnvType.T_C;
import static com.hartwig.hmftools.amber.purity.CanonicalSnvType.T_G;
import static com.hartwig.hmftools.common.amber.AmberBase.A;
import static com.hartwig.hmftools.common.amber.AmberBase.C;
import static com.hartwig.hmftools.common.amber.AmberBase.G;
import static com.hartwig.hmftools.common.amber.AmberBase.T;

import org.junit.Assert;
import org.junit.Test;

public class CanonicalSnvTypeTest
{
    @Test
    public void typeTest()
    {
        Assert.assertEquals(C_A, CanonicalSnvType.type(C, A));
        Assert.assertEquals(C_A, CanonicalSnvType.type(G, T));
        Assert.assertEquals(C_G, CanonicalSnvType.type(C, G));
        Assert.assertEquals(C_G, CanonicalSnvType.type(G, C));
        Assert.assertEquals(C_T, CanonicalSnvType.type(C, T));
        Assert.assertEquals(C_T, CanonicalSnvType.type(G, A));
        Assert.assertEquals(T_A, CanonicalSnvType.type(T, A));
        Assert.assertEquals(T_A, CanonicalSnvType.type(A, T));
        Assert.assertEquals(T_C, CanonicalSnvType.type(T, C));
        Assert.assertEquals(T_C, CanonicalSnvType.type(A, G));
        Assert.assertEquals(T_G, CanonicalSnvType.type(T, G));
        Assert.assertEquals(T_G, CanonicalSnvType.type(A, C));
    }
}
