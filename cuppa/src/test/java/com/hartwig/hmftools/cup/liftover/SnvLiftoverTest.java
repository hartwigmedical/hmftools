package com.hartwig.hmftools.cup.liftover;

import org.junit.Test;
import static org.junit.Assert.assertEquals;

import java.io.File;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

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
    public void canRunSnvLiftover()
    {
        TMP_DIR.mkdir();

        SnvLiftover liftover = new SnvLiftover(buildConfig());
        liftover.run();

        for(File file : TMP_DIR.listFiles()) file.delete();
        TMP_DIR.delete();
    }
}
