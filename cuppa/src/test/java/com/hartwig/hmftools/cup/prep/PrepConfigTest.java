package com.hartwig.hmftools.cup.prep;

import static com.hartwig.hmftools.cup.CuppaConfig.DNA_CATEGORIES;
import static com.hartwig.hmftools.cup.common.CupConstants.APP_NAME;

import static org.junit.Assert.assertEquals;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.junit.Test;

public class PrepConfigTest
{
    @Test
    public void canParseArgsMultiSample()
    {
        String sampleDataDir = Resources.getResource("pipeline_output/").getPath();
        String sampleIdFile = Resources.getResource("pipeline_output/sample_ids.csv").getPath();
        String refAltSjSites = Resources.getResource("alt_sj.selected_loci.minimal.tsv").getPath();
        String dummyValidPath = "/tmp/";

        String[] args = {
                // "-sample", sample,
                "-sample_id_file", sampleIdFile,
                "-categories", DNA_CATEGORIES,
                "-ref_genome_version", "V37",

                // Wild-cards are required to get subdirs correctly
                "-sample_data_dir", sampleDataDir + "*",
                "-linx_dir", sampleDataDir + "*/linx/",
                "-purple_dir", sampleDataDir + "*/purple/",
                //"-virus_dir", sampleDataDir + "*/virus_interpreter/",

                "-output_dir", dummyValidPath,
                "-isofox_dir", dummyValidPath,
                "-ref_alt_sj_sites", refAltSjSites,

                "-write_by_category"
        };

        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        PrepConfig.registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);
        PrepConfig prepConfig = new PrepConfig(configBuilder);

        assertEquals(2, prepConfig.SampleIds.size());
        assertEquals("V37", prepConfig.RefGenVersion.toString());
        assertEquals(true, prepConfig.WriteByCategory);
        assertEquals(refAltSjSites, prepConfig.AltSpliceJunctionSites);

        String selectedSampleId = "SKINMERKEL01T";
        assertEquals(
                sampleDataDir + selectedSampleId + "/linx/",
                prepConfig.getLinxDataDir(selectedSampleId)
        );

        assertEquals(
                sampleDataDir + selectedSampleId + "/purple/",
                prepConfig.getPurpleDataDir(selectedSampleId)
        );

        // When "-{tool_name}_dir" is not specified, the sample dir is used
        assertEquals(
                sampleDataDir + selectedSampleId,
                prepConfig.getVirusDataDir(selectedSampleId)
        );
    }
}
