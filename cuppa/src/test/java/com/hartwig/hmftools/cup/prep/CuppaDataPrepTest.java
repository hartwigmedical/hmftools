package com.hartwig.hmftools.cup;
import static com.hartwig.hmftools.cup.common.CupConstants.APP_NAME;

import static org.junit.Assert.assertTrue;

import java.io.File;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.cup.prep.CuppaDataPrep;
import com.hartwig.hmftools.cup.prep.PrepConfig;

import org.junit.Test;

public class CuppaDataPrepTest
{
    @Test
    public void canExtractFeaturesFromOneSample()
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        PrepConfig.registerConfig(configBuilder);

        String SAMPLE_DATA_DIR = Resources.getResource("pipeline_output/prostate_sample/").getPath();

        String[] args = {
                "-sample", "prostate_sample",
                "-categories", "DNA",
                "-ref_genome_version", "V37",
                "-sample_data_dir", SAMPLE_DATA_DIR,
                "-output_dir", "/Users/lnguyen/Desktop/test_output/" // TODO: Change to tmp path
        };

        configBuilder.checkAndParseCommandLine(args);

        CuppaDataPrep cuppaDataPrep = new CuppaDataPrep(configBuilder);
        cuppaDataPrep.runSingleSample();

        File outputFile = new File(cuppaDataPrep.mSampleDataWriter.getOutputPathSingleSample());
        assertTrue(outputFile.exists());
        outputFile.delete();
    }

    @Test
    public void canExtractFeaturesFromManySamples()
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        PrepConfig.registerConfig(configBuilder);

        String SAMPLE_DATA_DIR = Resources.getResource("pipeline_output/").getPath();

        String[] args = {
                "-sample", "prostate_sample",
                "-categories", "DNA",
                "-ref_genome_version", "V37",
                "-sample_data_dir", SAMPLE_DATA_DIR,
                "-output_dir", "/Users/lnguyen/Desktop/test_output/" // TODO: Change to tmp path
        };

        configBuilder.checkAndParseCommandLine(args);

        CuppaDataPrep cuppaDataPrep = new CuppaDataPrep(configBuilder);
        cuppaDataPrep.runMultiSample();

//        File outputFile = new File(cuppaDataPrep.mSampleDataWriter.getOutputPathSingleSample());
//        assertTrue(outputFile.exists());
//        outputFile.delete();
    }
}
