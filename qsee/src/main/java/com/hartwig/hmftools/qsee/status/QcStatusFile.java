package com.hartwig.hmftools.qsee.status;

import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_NAME;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_TYPE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_VALUE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SAMPLE_ID;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SAMPLE_TYPE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SOURCE_TOOL;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.MULTISAMPLE_SAMPLE_ID;
import static com.hartwig.hmftools.qsee.feature.FeatureMetadata.FIELD_QC_STATUS;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;
import com.hartwig.hmftools.qsee.QseeConfig;
import com.hartwig.hmftools.qsee.common.QseeFileCommon;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.NumberFormat;
import com.hartwig.hmftools.qsee.feature.SampleFeatures;

import org.jetbrains.annotations.Nullable;

public class QcStatusFile
{
    public static final String COL_QC_STATUS = FIELD_QC_STATUS;
    public static final String COL_FAIL_CONDITION = "FailCondition";

    public static String generateFilename(String basePath, String sampleId, @Nullable String outputId)
    {
        return QseeFileCommon.generateFilename(basePath, sampleId, "status", outputId, "tsv.gz");
    }

    public static String generateFilename(QseeConfig config)
    {
        String sampleId = config.isSinglePatient() ?
                config.getSampleIds(SampleType.TUMOR).get(0) :
                MULTISAMPLE_SAMPLE_ID;

        return generateFilename(config.OutputDir, sampleId, config.OutputId);
    }

    private record SampleFeatureRow(String sampleId, SampleType sampleType, Feature feature) {}

    public static void write(String filename, List<SampleFeatures> multiSampleFeatures)
    {
        List<SampleFeatureRow> rowsWithQcStatus = new ArrayList<>();
        for(SampleFeatures sampleFeatures : multiSampleFeatures)
        {
            for(Feature feature : sampleFeatures.features())
            {
                QcStatusType qcStatusType = feature.metadata().qcStatus().type();

                if(qcStatusType != QcStatusType.NONE)
                    rowsWithQcStatus.add(new SampleFeatureRow(sampleFeatures.sampleId(), sampleFeatures.sampleType(), feature));
            }
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

        DelimFileWriter.write(filename, columns, rowsWithQcStatus, (sampleFeatureRow, row) -> {
            row.set(COL_SAMPLE_ID, sampleFeatureRow.sampleId());
            row.set(COL_SAMPLE_TYPE, sampleFeatureRow.sampleType().name());

            FeatureKey featureKey = sampleFeatureRow.feature().key();
            row.set(COL_SOURCE_TOOL, featureKey.sourceTool().name());
            row.set(COL_FEATURE_TYPE, featureKey.type().name());
            row.set(COL_FEATURE_NAME, featureKey.name());
            row.set(COL_FEATURE_VALUE, sampleFeatureRow.feature().value());

            QcStatus qcStatus = sampleFeatureRow.feature().metadata().qcStatus();
            row.set(COL_QC_STATUS, qcStatus.type().name());
            row.set(COL_FAIL_CONDITION, qcStatus.displayString(NumberFormat.NUMBER));
        });
    }
}
