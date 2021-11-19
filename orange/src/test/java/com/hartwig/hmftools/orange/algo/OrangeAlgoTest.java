package com.hartwig.hmftools.orange.algo;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.orange.ImmutableOrangeConfig;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.OrangeConfigTestFactory;

import org.junit.Test;

public class OrangeAlgoTest {

    @Test
    public void canCreateReportFromTestDir() throws IOException {
        OrangeConfig config = OrangeConfigTestFactory.createTestOrangeConfig();
        OrangeAlgo algo = OrangeAlgo.fromConfig(config);

        assertNotNull(algo.run(config));
    }

    @Test
    public void canCreateReportWithoutTumorDoids() throws IOException {
        OrangeConfig config = ImmutableOrangeConfig.builder()
                .from(OrangeConfigTestFactory.createTestOrangeConfig())
                .primaryTumorDoids(Sets.newHashSet())
                .build();

        OrangeAlgo algo = OrangeAlgo.fromConfig(config);

        assertNotNull(algo.run(config));
    }
}