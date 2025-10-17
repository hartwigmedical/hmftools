package com.hartwig.hmftools.common.pipeline;

import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.DB_V6_0_FORMAT;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.OA_V2_0_FORMAT;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.OA_V2_2_FORMAT;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.PIP5_V6_0_FORMAT;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.PIPELINE_FORMAT_CFG;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.PIPELINE_FORMAT_FILE_CFG;

import static org.junit.Assert.assertEquals;

import java.io.File;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.junit.Test;

public class PipelineToolDirectoriesTest
{
    private static final String RESOURCE_DIR = Resources.getResource("pipeline").getPath();
    private static final String PARTIAL_TEST_CONFIG_FILE = RESOURCE_DIR + File.separator + "partialToolDirectoryConfig.tsv";
    private static final String COMPLETE_TEST_CONFIG_FILE = RESOURCE_DIR + File.separator + "completeToolDirectoryConfig.tsv";

    @Test
    public void defaultsToPipeline5Structure()
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        PipelineToolDirectories.addPipelineFormatOptions(configBuilder);

        PipelineToolDirectories victim = PipelineToolDirectories.resolveToolDirectories(
                configBuilder, PIPELINE_FORMAT_CFG, PIPELINE_FORMAT_FILE_CFG);
        assertEquals("pave", victim.paveGermlineDir());
        assertEqualDirectories(OA_V2_2_FORMAT, victim);
    }

    @Test
    public void canSpecifyPipeline5Structure()
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        PipelineToolDirectories.addPipelineFormatOptions(configBuilder);
        setPipelineConfig(configBuilder, PipelineOutputStructure.PIP5_V6_0.toString());

        PipelineToolDirectories victim = PipelineToolDirectories.resolveToolDirectories(
                configBuilder, PIPELINE_FORMAT_CFG, PIPELINE_FORMAT_FILE_CFG);
        assertEquals("pave_germline", victim.paveGermlineDir());
        assertEqualDirectories(PIP5_V6_0_FORMAT, victim);
    }

    private static void setPipelineConfig(final ConfigBuilder configBuilder, final String pipelineVersion)
    {
        configBuilder.checkAndParseCommandLine(new String[]{"-" + PIPELINE_FORMAT_CFG, pipelineVersion});
    }

    @Test
    public void canSpecifyOncoAnalyserStructure()
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        PipelineToolDirectories.addPipelineFormatOptions(configBuilder);
        setPipelineConfig(configBuilder, PipelineOutputStructure.OA_V2_2.toString());

        PipelineToolDirectories victim = PipelineToolDirectories.resolveToolDirectories(
                configBuilder, PIPELINE_FORMAT_CFG, PIPELINE_FORMAT_FILE_CFG);
        assertEquals("pave", victim.paveGermlineDir());
        assertEqualDirectories(OA_V2_2_FORMAT, victim);
    }

    @Test
    public void canSpecifyDatabaseStructure()
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        PipelineToolDirectories.addPipelineFormatOptions(configBuilder);
        setPipelineConfig(configBuilder, PipelineOutputStructure.DB_V6_0.toString());

        PipelineToolDirectories victim = PipelineToolDirectories.resolveToolDirectories(configBuilder, PIPELINE_FORMAT_CFG, PIPELINE_FORMAT_FILE_CFG);
        assertEquals("pave/germline", victim.paveGermlineDir());
        assertEqualDirectories(DB_V6_0_FORMAT, victim);
    }

    @Test
    public void canLoadFromPartialFile()
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        PipelineToolDirectories.addPipelineFormatOptions(configBuilder);
        configBuilder.checkAndParseCommandLine(new String[]{"-" + PIPELINE_FORMAT_FILE_CFG, PARTIAL_TEST_CONFIG_FILE });

        PipelineToolDirectories victim = PipelineToolDirectories.resolveToolDirectories(configBuilder, PIPELINE_FORMAT_CFG, PIPELINE_FORMAT_FILE_CFG);
        assertEquals("amberTest", victim.amberDir());
        assertEquals("purple/*/test", victim.purpleDir());
        assertEquals("sage/$", victim.sageGermlineDir());
        assertEquals("", victim.sageSomaticDir());
    }

    @Test
    public void canLoadFromCompleteFile()
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        PipelineToolDirectories.addPipelineFormatOptions(configBuilder);
        configBuilder.checkAndParseCommandLine(new String[]{"-" + PIPELINE_FORMAT_FILE_CFG, COMPLETE_TEST_CONFIG_FILE });

        PipelineToolDirectories victim =
                PipelineToolDirectories.resolveToolDirectories(configBuilder, PIPELINE_FORMAT_CFG, PIPELINE_FORMAT_FILE_CFG);
        assertEquals(OA_V2_0_FORMAT, victim);
    }

    @Test
    public void preferFileOverOtherArguments()
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        PipelineToolDirectories.addPipelineFormatOptions(configBuilder);
        String[] arguments = { "-" + PIPELINE_FORMAT_CFG, "DB_V6_0", "-" + PIPELINE_FORMAT_FILE_CFG, PARTIAL_TEST_CONFIG_FILE };
        configBuilder.checkAndParseCommandLine(arguments);

        PipelineToolDirectories victim =
                PipelineToolDirectories.resolveToolDirectories(configBuilder, PIPELINE_FORMAT_CFG, PIPELINE_FORMAT_FILE_CFG);
        assertEquals("amberTest", victim.amberDir());
        assertEquals("purple/*/test", victim.purpleDir());
        assertEquals("sage/$", victim.sageGermlineDir());
        assertEquals("", victim.sageSomaticDir());
    }

    @Test
    public void canGetWildcardPaths()
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        PipelineToolDirectories.addPipelineFormatOptions(configBuilder);

        setPipelineConfig(configBuilder, PipelineOutputStructure.PIP5_V6_0.toString());
        PipelineToolDirectories victim = PipelineToolDirectories.resolveToolDirectories(
                configBuilder, PIPELINE_FORMAT_CFG, PIPELINE_FORMAT_FILE_CFG);
        assertEquals("*/bam_metrics", victim.tumorMetricsDir());
        assertEquals("$/bam_metrics", victim.germlineMetricsDir());
        assertEquals("purple", victim.purpleDir());
        assertEqualDirectories(PIP5_V6_0_FORMAT, victim);
    }

    @Test
    public void canResolveTumorSampleId()
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        PipelineToolDirectories.addPipelineFormatOptions(configBuilder);
        setPipelineConfig(configBuilder, PipelineOutputStructure.PIP5_V6_0.toString());

        String tumorSampleId = "TUMOR";

        PipelineToolDirectories victim = PipelineToolDirectories.resolveToolDirectories(
                configBuilder, PIPELINE_FORMAT_CFG, PIPELINE_FORMAT_FILE_CFG, tumorSampleId);
        assertEquals("TUMOR/bam_metrics", victim.tumorMetricsDir());
        assertEquals("$/bam_metrics", victim.germlineMetricsDir());
        assertEquals("purple", victim.purpleDir());
    }

    @Test
    public void canResolveBothSampleIds()
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        PipelineToolDirectories.addPipelineFormatOptions(configBuilder);
        setPipelineConfig(configBuilder, PipelineOutputStructure.PIP5_V6_0.toString());

        String tumorSampleId = "TUMOR";
        String normalSampleId = "NORMAL";

        PipelineToolDirectories victim = PipelineToolDirectories.resolveToolDirectories(
                configBuilder, PIPELINE_FORMAT_CFG, PIPELINE_FORMAT_FILE_CFG, tumorSampleId, normalSampleId);
        assertEquals("TUMOR/bam_metrics", victim.tumorMetricsDir());
        assertEquals("NORMAL/bam_metrics", victim.germlineMetricsDir());
        assertEquals("purple", victim.purpleDir());
    }

    private void assertEqualDirectories(final PipelineToolDirectories expected, final PipelineToolDirectories victim)
    {
        assertEquals(expected.amberDir(), victim.amberDir());
        assertEquals(expected.chordDir(), victim.chordDir());
        assertEquals(expected.ciderDir(), victim.ciderDir());
        assertEquals(expected.cobaltDir(), victim.cobaltDir());
        assertEquals(expected.cuppaDir(), victim.cuppaDir());
        assertEquals(expected.esveeDir(), victim.esveeDir());
        assertEquals(expected.germlineFlagstatDir(), victim.germlineFlagstatDir());
        assertEquals(expected.germlineMetricsDir(), victim.germlineMetricsDir());
        assertEquals(expected.isofoxDir(), victim.isofoxDir());
        assertEquals(expected.lilacDir(), victim.lilacDir());
        assertEquals(expected.linxGermlineDir(), victim.linxGermlineDir());
        assertEquals(expected.linxSomaticDir(), victim.linxSomaticDir());
        assertEquals(expected.orangeDir(), victim.orangeDir());
        assertEquals(expected.paveGermlineDir(), victim.paveGermlineDir());
        assertEquals(expected.paveSomaticDir(), victim.paveSomaticDir());
        assertEquals(expected.peachDir(), victim.peachDir());
        assertEquals(expected.purpleDir(), victim.purpleDir());
        assertEquals(expected.sageGermlineDir(), victim.sageGermlineDir());
        assertEquals(expected.sageSomaticDir(), victim.sageSomaticDir());
        assertEquals(expected.sigsDir(), victim.sigsDir());
        assertEquals(expected.snpGenotypeDir(), victim.snpGenotypeDir());
        assertEquals(expected.tealDir(), victim.tealDir());
        assertEquals(expected.tumorFlagstatDir(), victim.tumorFlagstatDir());
        assertEquals(expected.tumorMetricsDir(), victim.tumorMetricsDir());
        assertEquals(expected.virusBreakendDir(), victim.virusBreakendDir());
        assertEquals(expected.virusInterpreterDir(), victim.virusInterpreterDir());
    }
}
