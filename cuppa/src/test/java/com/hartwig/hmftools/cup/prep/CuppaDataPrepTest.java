package com.hartwig.hmftools.cup.prep;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.CuppaConfig.DNA_CATEGORIES;
import static com.hartwig.hmftools.cup.common.CupConstants.APP_NAME;

import static org.junit.Assert.assertTrue;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.junit.Test;

public class CuppaDataPrepTest
{

    @Test
    public void canExtractFeaturesSingleSample()
    {
        String sampleDataDir = Resources.getResource("pipeline_output/SKINMERKEL01T/").getPath();

        String[] args = {
                // "-sample", "PROSTATE01T",
                "-sample", "SKINMERKEL01T",
                "-categories", "DNA",
                "-ref_genome_version", "V37",
                "-sample_data_dir", sampleDataDir,
                "-output_dir", "/Users/lnguyen/Desktop/test_output/" // TODO: Change to tmp path
        };
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        PrepConfig.registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        CuppaDataPrep cuppaDataPrep = new CuppaDataPrep(configBuilder);
        List<DataItem> dataItems = cuppaDataPrep.getDataSingleSample();
        CuppaDataPrep.writeDataSingleSample(dataItems, cuppaDataPrep.mConfig.OutputDir + "features.tsv.gz");

//        File outputFile = new File(cuppaDataPrep.mSampleDataWriter.getOutputPathSingleSample());
//        assertTrue(outputFile.exists());
//        outputFile.delete();
    }

    @Test
    public void canGetDataItemMatrix() throws IOException
    {
        String sampleDataDir = Resources.getResource("pipeline_output/").getPath();
        String sampleIdFile = Resources.getResource("pipeline_output/sample_ids.csv").getPath();
        // String refAltSjSites = Resources.getResource("alt_sj.selected_loci.minimal.tsv").getPath();

        String[] args = {
                // "-sample", sample,
                "-sample_id_file", sampleIdFile,
                "-categories", DNA_CATEGORIES,
                "-ref_genome_version", "V37",

                "-sample_data_dir", sampleDataDir + "*", // Wild-cards are required to get subdirs correctly
                "-output_dir",  "/Users/lnguyen/Desktop/test_output/",
                // "-isofox_dir", dummyValidPath,
                // "-ref_alt_sj_sites", refAltSjSites,

                // "-write_by_category"
        };

        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        PrepConfig.registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        CuppaDataPrep cuppaDataPrep = new CuppaDataPrep(configBuilder);
//        cuppaDataPrep.run();

//        System.out.println(cuppaDataPrep.mDataPreparers.get(3).toString());
        DataItemMatrix dataItemMatrix = cuppaDataPrep.getDataOneCategoryMultiSample(
                cuppaDataPrep.mDataPreparers.get(3)
        );

        CuppaDataPrep.writeDataMultiSample(dataItemMatrix, "/Users/lnguyen/Desktop/test_output/features.tsv.gz", false);
    }
}
