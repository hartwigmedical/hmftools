package com.hartwig.hmftools.orange.algo;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.orange.ImmutableOrangeConfig;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.TestOrangeConfigFactory;

import org.junit.Test;

public class OrangeAlgoTest {

    @Test
    public void canRunReportFromTestDirDNATumorOnly() throws IOException {
        OrangeConfig config = TestOrangeConfigFactory.createDNAConfigTumorOnly();
        OrangeAlgo algo = OrangeAlgo.fromConfig(config);

        assertNotNull(algo.run(config));
    }

    @Test
    public void canRunReportFromTestDirDNA() throws IOException {
        OrangeConfig config = TestOrangeConfigFactory.createDNAConfigTumorNormal();
        OrangeAlgo algo = OrangeAlgo.fromConfig(config);

        assertNotNull(algo.run(config));
    }

    @Test
    public void canRunReportFromTestDirDNARNA() throws IOException {
        OrangeConfig config = TestOrangeConfigFactory.createDNARNAConfig();
        OrangeAlgo algo = OrangeAlgo.fromConfig(config);

        assertNotNull(algo.run(config));
    }

    @Test
    public void canCreateReportWithoutTumorDoids() throws IOException {
        OrangeConfig config = ImmutableOrangeConfig.builder()
                .from(TestOrangeConfigFactory.createDNAConfigTumorNormal())
                .primaryTumorDoids(Sets.newHashSet())
                .build();

        OrangeAlgo algo = OrangeAlgo.fromConfig(config);

        assertNotNull(algo.run(config));
    }
}