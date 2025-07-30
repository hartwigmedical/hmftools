package com.hartwig.hmftools.common.utils.r;

import java.io.IOException;

import com.hartwig.hmftools.common.utils.RExecutor;

import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class RExecutorTest
{
    @Test
    public void testR() throws IOException, InterruptedException
    {
        RExecutor.executeFromClasspath("r/dummyR.R", "0");
    }
}
