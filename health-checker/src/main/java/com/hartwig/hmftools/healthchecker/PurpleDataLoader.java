package com.hartwig.hmftools.healthchecker;

import static com.hartwig.hmftools.healthchecker.HealthChecksApplication.HC_LOGGER;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.purple.PurpleQCFile;

public class PurpleDataLoader
{
    public static List<QCValue> loadPurpleValues(final String sampleId, final String purpleDir)
    {
        String qcFilename = PurpleQCFile.generateFilename(purpleDir, sampleId);

        try
        {
            PurpleQC qcData = PurpleQCFile.read(qcFilename);

            return Lists.newArrayList(
                    new QCValue(QCValueType.PURPLE_QC_STATUS, qcData.toString()),
                    new QCValue(QCValueType.PURPLE_CONTAMINATION, String.valueOf(qcData.contamination())));
        }
        catch(IOException e)
        {
            HC_LOGGER.error("failed to read purple QC from ({}): {}", qcFilename, e.toString());
            return null;
        }
    }
}
