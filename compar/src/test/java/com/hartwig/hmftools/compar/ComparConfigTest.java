package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.compar.ComparConfig.TRUTHSET_FILES;
import static com.hartwig.hmftools.compar.common.ComparConstants.FLD_CATEGORY;
import static com.hartwig.hmftools.compar.common.ComparConstants.FLD_ITEM_KEY;
import static com.hartwig.hmftools.compar.common.SourceType.NEW;
import static com.hartwig.hmftools.compar.common.SourceType.OLD;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.SourceType;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

public class ComparConfigTest
{
    @Rule
    public TemporaryFolder TempDir = new TemporaryFolder();

    @Test
    public void testSampleDirColumnsOverrideDefaultSampleDir() throws IOException
    {
        ComparConfig config = createConfig(List.of(
                "SampleId,OldSampleDir,NewSampleDir",
                "SAMPLE_1,/run_a/SAMPLE_1,/run_b/SAMPLE_1/",
                "SAMPLE_2,,"));

        assertTrue(config.isValid());

        assertEquals("/run_a/SAMPLE_1/purple/", purpleDir(config, OLD, "SAMPLE_1"));
        assertEquals("/run_b/SAMPLE_1/purple/", purpleDir(config, NEW, "SAMPLE_1"));

        // empty values fall back to the default sample directory
        assertEquals("/default_old/purple/", purpleDir(config, OLD, "SAMPLE_2"));
        assertEquals("/default_new/purple/", purpleDir(config, NEW, "SAMPLE_2"));
    }

    @Test
    public void testSampleDirColumnsAreIndependentPerSource() throws IOException
    {
        ComparConfig config = createConfig(List.of(
                "SampleId,OldSampleDir",
                "SAMPLE_1,/run_a/SAMPLE_1/"));

        assertTrue(config.isValid());

        assertEquals("/run_a/SAMPLE_1/purple/", purpleDir(config, OLD, "SAMPLE_1"));
        assertEquals("/default_new/purple/", purpleDir(config, NEW, "SAMPLE_1"));
    }

    @Test
    public void testSampleDirAppliesToSourceSampleId() throws IOException
    {
        ComparConfig config = createConfig(List.of(
                "SampleId,OldSampleId,NewSampleId,OldSampleDir",
                "SAMPLE_1,SAMPLE_1_OLD,SAMPLE_1_NEW,/run_a/SAMPLE_1_OLD/"));

        assertTrue(config.isValid());

        assertEquals("/run_a/SAMPLE_1_OLD/purple/", purpleDir(config, OLD, "SAMPLE_1"));
        assertEquals("/default_new/purple/", purpleDir(config, NEW, "SAMPLE_1"));
    }

    @Test
    public void testSampleDirWithWildcardToolDir() throws IOException
    {
        ConfigBuilder configBuilder = createConfigBuilder(List.of(
                "SampleId,OldSampleDir",
                "SAMPLE_1,/run_a/"));

        configBuilder.setValue(PURPLE_DIR_CFG + "_old", "*/purple/");

        ComparConfig config = new ComparConfig(configBuilder);

        assertTrue(config.isValid());
        assertEquals("/run_a/SAMPLE_1/purple/", purpleDir(config, OLD, "SAMPLE_1"));
    }

    @Test
    public void testSampleDirForTruthsetSourceIsInvalid() throws IOException
    {
        File truthsetFile = TempDir.newFile("truthset.tsv");
        Files.write(truthsetFile.toPath(), List.of(FLD_CATEGORY + "\t" + FLD_ITEM_KEY));

        ConfigBuilder configBuilder = createConfigBuilder(List.of(
                "SampleId,OldSampleDir",
                "SAMPLE_1,/run_a/SAMPLE_1/"));

        configBuilder.setValue(TRUTHSET_FILES + "_old", truthsetFile.getAbsolutePath());

        ComparConfig config = new ComparConfig(configBuilder);

        assertFalse(config.isValid());
    }

    private ComparConfig createConfig(final List<String> sampleIdLines) throws IOException
    {
        return new ComparConfig(createConfigBuilder(sampleIdLines));
    }

    private ConfigBuilder createConfigBuilder(final List<String> sampleIdLines) throws IOException
    {
        File sampleIdFile = TempDir.newFile("sample_ids.csv");
        Files.write(sampleIdFile.toPath(), sampleIdLines);

        ConfigBuilder configBuilder = new ConfigBuilder();
        ComparConfig.addConfig(configBuilder);

        configBuilder.setValue(SAMPLE_ID_FILE, sampleIdFile.getAbsolutePath());
        configBuilder.setValue("sample_dir_old", "/default_old");
        configBuilder.setValue("sample_dir_new", "/default_new");

        // explicit tool directories keep the expected paths independent of the default pipeline format
        configBuilder.setValue(PURPLE_DIR_CFG + "_old", "purple/");
        configBuilder.setValue(PURPLE_DIR_CFG + "_new", "purple/");

        return configBuilder;
    }

    private static String purpleDir(final ComparConfig config, final SourceType sourceType, final String sampleId)
    {
        PipelineSourcePaths sourcePaths = PipelineSourcePaths.sampleInstance(
                config.getSourceData(sourceType).PipelinePaths,
                config.sourceSampleId(sourceType, sampleId),
                config.sourceReferenceId(sourceType, sampleId));

        return sourcePaths.Purple;
    }
}
