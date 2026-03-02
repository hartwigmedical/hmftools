package com.hartwig.hmftools.qsee.prep;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_NAME;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_TYPE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_VALUE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SAMPLE_ID;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SAMPLE_TYPE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SOURCE_TOOL;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.qsee.common.QseeFileCommon;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureKey;

import org.jetbrains.annotations.Nullable;

public class VisDataFile
{
    private static final String COL_PLOT_METADATA = "PlotMetadata";

    public static String generateFilename(String basePath, String sampleId, @Nullable String outputId)
    {
        return QseeFileCommon.generateFilename(basePath, sampleId, "vis.data", outputId, "tsv.gz");
    }

    public static String generateFilename(QseePrepConfig config)
    {
        String sampleId = config.isSinglePatient() ?
                config.getSampleIds(SampleType.TUMOR).get(0) :
                "multisample";

        return generateFilename(config.OutputDir, sampleId, config.OutputId);
    }

    public static void write(String outputFile, List<VisSampleData> visSampleDataEntries)
    {
        try(BufferedWriter writer = createBufferedWriter(outputFile))
        {
            QC_LOGGER.info("Writing sample vis data to: {}", outputFile);

            StringJoiner header = new StringJoiner(TSV_DELIM);

            header.add(COL_SAMPLE_ID);
            header.add(COL_SAMPLE_TYPE);
            header.add(COL_SOURCE_TOOL);
            header.add(COL_FEATURE_TYPE);
            header.add(COL_FEATURE_NAME);
            header.add(COL_FEATURE_VALUE);
            header.add(COL_PLOT_METADATA);

            writer.write(header.toString());
            writer.newLine();

            for(VisSampleData entry : visSampleDataEntries)
            {
                Feature feature = entry.feature();
                FeatureKey featureKey = feature.key();

                StringJoiner line = new StringJoiner(TSV_DELIM);
                line.add(entry.sampleId());
                line.add(entry.sampleType().name());
                line.add(featureKey.sourceTool().name());
                line.add(featureKey.type().name());
                line.add(featureKey.name());

                String featureValue = QseeFileCommon.DECIMAL_FORMAT.format(feature.value());
                line.add(featureValue);

                line.add(feature.metadata().displayString());

                writer.write(line.toString());
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            QC_LOGGER.error("Failed to write to file: {}", outputFile, e);
            System.exit(1);
        }
    }
}
