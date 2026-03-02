package com.hartwig.hmftools.qsee.status;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.qsee.common.QseeConstants.APP_NAME;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class DefaultThresholdsWriter
{
    private static final String NO_OVERRIDABLE_ONLY = "no_overridable_only";
    private static final String NO_OVERRIDABLE_ONLY_DESC = "Also write thresholds that are not overridable";

    private static String generateFilename(String basePath)
    {
        return checkAddDirSeparator(basePath) + "qsee.thresholds.default.tsv";
    }

    private static void writeDefaults(String filename, boolean overridableOnly)
    {
        ThresholdRegistry defaultThresholds = ThresholdRegistry.createDefault();

        List<QcThreshold> thresholdsToWrite = new ArrayList<>(defaultThresholds.getAll());

        thresholdsToWrite.sort(Comparator.comparing(threshold -> threshold.key().sampleType()));

        if(overridableOnly)
        {
            QC_LOGGER.info("Skipping thresholds that cannot be overridden");

            thresholdsToWrite = thresholdsToWrite.stream()
                    .filter(threshold -> !threshold.determinedElsewhere())
                    .toList();
        }

        ThresholdOverridesFile.write(filename, thresholdsToWrite);
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addPath(OUTPUT_DIR, false, OUTPUT_DIR);
        configBuilder.addFlag(NO_OVERRIDABLE_ONLY, NO_OVERRIDABLE_ONLY_DESC);
        configBuilder.checkAndParseCommandLine(args);

        String outputDir = configBuilder.getValue(OUTPUT_DIR, "");
        String outputFile = generateFilename(outputDir);

        boolean overridableOnly = !configBuilder.hasFlag(NO_OVERRIDABLE_ONLY);

        QC_LOGGER.info("Writing default thresholds to: {}", outputFile);
        writeDefaults(outputFile, overridableOnly);
    }
}
