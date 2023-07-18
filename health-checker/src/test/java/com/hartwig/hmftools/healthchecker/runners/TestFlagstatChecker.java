package com.hartwig.hmftools.healthchecker.runners;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.flagstat.Flagstat;
import com.hartwig.hmftools.common.flagstat.FlagstatFile;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;

import org.jetbrains.annotations.Nullable;

public class TestFlagstatChecker implements HealthChecker
{
    private final String mRefFlagstatFile;
    @Nullable
    private final String mTumorFlagstatFile;

    public TestFlagstatChecker(final String refFlagstatFile, @Nullable final String tumorFlagstatFile)
    {
        mRefFlagstatFile = refFlagstatFile;
        mTumorFlagstatFile = tumorFlagstatFile;
    }

    @Override
    public List<QCValue> run() throws IOException
    {
        List<QCValue> qcValues = Lists.newArrayList();

        Flagstat refFlagstat = FlagstatFile.read(mRefFlagstatFile);
        qcValues.add(new QCValue(QCValueType.REF_PROPORTION_MAPPED, String.valueOf(refFlagstat.mappedProportion())));
        qcValues.add(new QCValue(QCValueType.REF_PROPORTION_DUPLICATE, String.valueOf(refFlagstat.duplicateProportion())));

        if(mTumorFlagstatFile != null)
        {
            Flagstat tumFlagstat = FlagstatFile.read(mTumorFlagstatFile);
            qcValues.add(new QCValue(QCValueType.TUM_PROPORTION_MAPPED, String.valueOf(tumFlagstat.mappedProportion())));
            qcValues.add(new QCValue(QCValueType.TUM_PROPORTION_DUPLICATE, String.valueOf(tumFlagstat.duplicateProportion())));
        }

        return qcValues;
    }
}
