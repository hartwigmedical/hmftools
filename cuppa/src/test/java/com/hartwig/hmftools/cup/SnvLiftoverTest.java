package com.hartwig.hmftools.cup;

import org.apache.commons.io.FileUtils;
import org.junit.Test;
import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.cup.liftover.LiftoverConfig;
import com.hartwig.hmftools.cup.liftover.SnvLiftover;
import com.hartwig.hmftools.cup.liftover.VcfPositionConverter;

public class SnvLiftoverTest
{
    private static final String SAMPLE_ID = "TUMOR_SAMPLE";
    private static final String SAMPLE_VCF_DIR = Resources.getResource("pipeline_output/" + SAMPLE_ID).getPath();
    private static final String VCF_FILE = SAMPLE_VCF_DIR + SAMPLE_ID + ".purple.somatic.vcf.gz";

    File TMP_DIR = new File(System.getProperty("java.io.tmpdir") + "/SnvLiftover/");

    private ConfigBuilder buildConfig()
    {
            String[] args = {
                    "-sample", SAMPLE_ID,
                    "-sample_vcf_dir", SAMPLE_VCF_DIR,
                    "-output_dir", TMP_DIR.toString(),
                    "-log_debug",
                    "-threads", "1"
            };

            ConfigBuilder config = new ConfigBuilder();
            LiftoverConfig.addOptions(config);
            config.checkAndParseCommandLine(args);

            return config;
    }

    @Test
    public void canLiftOverPosition()
    {
        LiftoverConfig config = new LiftoverConfig(buildConfig());
        VcfPositionConverter converter = new VcfPositionConverter(SAMPLE_ID, VCF_FILE, config);

        assertEquals(761264, converter.convertPosition("1", 696644));
        assertEquals(55181378, converter.convertPosition("7", 55249071));
    }

    @Test
    public void canRunSnvLiftover() throws IOException
    {
        TMP_DIR.mkdir();

        SnvLiftover liftover = new SnvLiftover(buildConfig());
        liftover.run();

        FileUtils.deleteDirectory(TMP_DIR);
    }
}
