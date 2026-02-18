package com.hartwig.hmftools.orange.algo;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.TestOrangeConfigFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class OrangeAlgoTest
{
    @Test
    public void canRunReportFromTestDirOnMinimalConfig() throws Exception
    {
        OrangeConfig config = TestOrangeConfigFactory.createMinimalConfig();
        OrangeAlgo algo = createOrangeAlgo(config);

        assertNotNull(algo.run(config));
    }

    @Test
    public void canRunReportFromTestDirTargeted() throws Exception
    {
        OrangeConfig config = TestOrangeConfigFactory.createTargetedConfig();
        OrangeAlgo algo = createOrangeAlgo(config);

        assertNotNull(algo.run(config));
    }

    @Test
    public void canRunReportFromTestDirWGSTumorOnly() throws Exception
    {
        OrangeConfig config = TestOrangeConfigFactory.createWGSConfigTumorOnly();
        OrangeAlgo algo = createOrangeAlgo(config);

        assertNotNull(algo.run(config));
    }

    @Test
    public void canRunReportFromTestDirWGSTumorNormal() throws Exception
    {
        OrangeConfig config = TestOrangeConfigFactory.createWGSConfigTumorNormal();
        OrangeAlgo algo = createOrangeAlgo(config);

        assertNotNull(algo.run(config));
    }

    @Test
    public void canRunReportFromTestDirWGTSTumorNormal() throws Exception
    {
        OrangeConfig config = TestOrangeConfigFactory.createWGTSConfigTumorNormal();
        OrangeAlgo algo = createOrangeAlgo(config);

        assertNotNull(algo.run(config));
    }

    @NotNull
    private static OrangeAlgo createOrangeAlgo(@NotNull OrangeConfig config) throws IOException
    {
        OrangeAlgo algo = OrangeAlgo.fromConfig(config);
        return algo;
    }
}