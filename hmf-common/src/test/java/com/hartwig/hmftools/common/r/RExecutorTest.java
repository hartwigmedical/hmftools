package com.hartwig.hmftools.common.r;

import java.io.IOException;

import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class RExecutorTest {

    @Test
    public void testR() throws IOException, InterruptedException {
        RExecutor.executeFromClasspath("r/dummyR.R", "0");
    }
}
