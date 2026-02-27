package com.hartwig.hmftools.qsee.status;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_NAME;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_TYPE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_VALUE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SAMPLE_ID;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SAMPLE_TYPE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SOURCE_TOOL;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.QSEE_FILE_ID;
import static com.hartwig.hmftools.qsee.feature.FeatureMetadata.FIELD_QC_STATUS;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.NumberFormat;
import com.hartwig.hmftools.qsee.prep.QseePrepConfig;
import com.hartwig.hmftools.qsee.prep.VisSampleData;

import org.jetbrains.annotations.Nullable;

public class QcSummaryFile
{
    public static final String COL_QC_STATUS = FIELD_QC_STATUS;
    public static final String COL_FAIL_CONDITION = "FailCondition";

    public static String generateFilename(String basePath, String sampleId, @Nullable String outputId)
    {
        String filename = checkAddDirSeparator(basePath) + sampleId + "." + QSEE_FILE_ID + ".summary";

        if(outputId != null)
            filename += "." + outputId;

        filename += ".tsv.gz";

        return filename;
    }

    public static String generateFilename(QseePrepConfig config)
    {
        String sampleId = config.isSinglePatient() ?
                config.getSampleIds(SampleType.TUMOR).get(0) :
                "multisample";

        return generateFilename(config.OutputDir, sampleId, config.OutputId);
    }

    public static void write(String filename, List<VisSampleData> visData)
    {
        List<VisSampleData> visDataWithQcStatus = new ArrayList<>();
        for(VisSampleData visSampleData : visData)
        {
            QcStatusType qcStatusType = visSampleData.feature().metadata().qcStatus().type();

            if(qcStatusType != QcStatusType.NONE)
                visDataWithQcStatus.add(visSampleData);
        }

        List<String> columns = List.of(
                COL_SAMPLE_ID,
                COL_SAMPLE_TYPE,
                COL_SOURCE_TOOL,
                COL_FEATURE_TYPE,
                COL_FEATURE_NAME,
                COL_FEATURE_VALUE,
                COL_QC_STATUS,
                COL_FAIL_CONDITION
        );

        DelimFileWriter.write(filename, columns, visDataWithQcStatus, (visSampleData, row) -> {
            row.set(COL_SAMPLE_ID, visSampleData.sampleId());
            row.set(COL_SAMPLE_TYPE, visSampleData.sampleType().name());

            FeatureKey featureKey = visSampleData.feature().key();
            row.set(COL_SOURCE_TOOL, featureKey.sourceTool().name());
            row.set(COL_FEATURE_TYPE, featureKey.type().name());
            row.set(COL_FEATURE_NAME, featureKey.name());
            row.set(COL_FEATURE_VALUE, visSampleData.feature().value());

            QcStatus qcStatus = visSampleData.feature().metadata().qcStatus();
            row.set(COL_QC_STATUS, qcStatus.type().name());
            row.set(COL_FAIL_CONDITION, qcStatus.displayString(NumberFormat.NUMBER));
        });
    }
}
