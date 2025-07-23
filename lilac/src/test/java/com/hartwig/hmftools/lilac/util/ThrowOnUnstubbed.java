package com.hartwig.hmftools.lilac.util;

import java.util.Arrays;

import org.mockito.invocation.InvocationOnMock;
import org.mockito.stubbing.Answer;

public class ThrowOnUnstubbed implements Answer<Object>
{
    @Override
    public Object answer(final InvocationOnMock invocation)
    {
        final String errMsg = "Unstubbed method call to: "
                + invocation.getMethod().getName()
                + ", with args: "
                + Arrays.toString(invocation.getArguments());
        throw new RuntimeException(errMsg);
    }
}
