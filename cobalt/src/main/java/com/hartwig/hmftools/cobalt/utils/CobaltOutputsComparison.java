package com.hartwig.hmftools.cobalt.utils;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltConstants.APP_NAME;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class CobaltOutputsComparison
{
    public static void main(@NotNull final String[] args) throws IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addPath(ORIGINAL_VALUES_DIR, true, "Original values directory");
        configBuilder.addPath(COMPARISON_VALUES_DIR, true, "Comparison values directory");
        configBuilder.addFlag(SAMPLE, "Sample name");

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        CobaltOutputsComparison comparer = new CobaltOutputsComparison(configBuilder);
        comparer.run();
    }

    public static final String ORIGINAL_VALUES_DIR = "original_values_dir";
    public static final String COMPARISON_VALUES_DIR = "comparison_values_dir";
    private final String mSampleId;
    private final File mOriginalValuesDirectory;
    private final File mComparisonValuesDirectory;
    private final File mOutputDirectory;

    public CobaltOutputsComparison(final ConfigBuilder configBuilder)
    {
        mSampleId = configBuilder.getValue(SAMPLE);
        mOriginalValuesDirectory = configBuilder.getFile(ORIGINAL_VALUES_DIR);
        mComparisonValuesDirectory = configBuilder.getFile(COMPARISON_VALUES_DIR);
        mOutputDirectory = configBuilder.getFile(OUTPUT_DIR);
    }

    public void run() throws IOException
    {
        String originalRatiosFileName = CobaltRatioFile.generateFilename(mOriginalValuesDirectory.getAbsolutePath(), mSampleId);
        List<RawCobaltRatio> originalRatios = new RawCobaltRatioFile(originalRatiosFileName).read();
        String comparisonRatiosFileName = CobaltRatioFile.generateFilename(mComparisonValuesDirectory.getAbsolutePath(), mSampleId);
        List<RawCobaltRatio> comparisonRatios = new RawCobaltRatioFile(comparisonRatiosFileName).read();

        int length = originalRatios.size();
        if(length != comparisonRatios.size())
        {
            CB_LOGGER.warn("Cobalt ratios files have different lengths. Original size: {}. Comparison size: {}", originalRatios.size(), comparisonRatios.size());
            return;
        }
        double epsilon = 0.01;
        List<RawCobaltRatio> differences = new ArrayList<>();
        for(int i = 0; i < length; i++)
        {
            RawCobaltRatio originalRatio = originalRatios.get(i);
            RawCobaltRatio comparisonRatio = comparisonRatios.get(i);
            if(originalRatio.matchesPosition(comparisonRatio))
            {
                RawCobaltRatio difference = originalRatio.differences(comparisonRatio, epsilon);
                if(difference != null)
                {
                    differences.add(difference);
                    differences.add(originalRatio);
                    differences.add(comparisonRatio);
                }
            }
            else
            {
                CB_LOGGER.warn("Original and comparison ratio files have mis-matched entries at index {}.", i);
            }
        }
        CB_LOGGER.info("Number of differences in ratios: {}", differences.size() / 3);
        String differencesFileName = checkAddDirSeparator(mOutputDirectory.getAbsolutePath()) + mSampleId + ".diff.tsv.gz";
        RawCobaltRatioFile.write(differencesFileName, differences);
    }
}
