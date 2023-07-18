package com.hartwig.hmftools.healthchecker.runners;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.purple.PurpleQCFile;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;

public class PurpleChecker implements HealthChecker
{
    private final String mTumorSample;
    private final String mPurpleDir;

    public PurpleChecker(final String tumorSample, final String purpleDirectory)
    {
        mTumorSample = tumorSample;
        mPurpleDir = purpleDirectory;
    }

    @Override
    public List<QCValue> run() throws IOException
    {
        String path = PurpleQCFile.generateFilename(mPurpleDir, mTumorSample);
        PurpleQC qc = PurpleQCFile.read(path);

        return Lists.newArrayList(
                new QCValue(QCValueType.PURPLE_QC_STATUS, qc.toString()),
                new QCValue(QCValueType.PURPLE_CONTAMINATION, String.valueOf(qc.contamination())));
    }
}
