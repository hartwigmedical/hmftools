package com.hartwig.hmftools.orange.algo;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.OrangeTestFactory;

import org.junit.Test;

public class OrangeAlgoTest {

    @Test
    public void canCreateReportFromTestDir() throws IOException {
        OrangeConfig config = OrangeTestFactory.createTestOrangeConfig();
        OrangeAlgo algo = OrangeAlgo.fromConfig(config);

        assertNotNull(algo.run(config));
    }
}