package com.hartwig.hmftools.healthchecker.runners;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.metrics.WGSMetrics;
import com.hartwig.hmftools.common.metrics.WGSMetricsFile;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;

public class ReferenceMetricsChecker extends FileBasedHealthChecker<WGSMetrics>
{
    public ReferenceMetricsChecker(final String metricsFile) throws IOException
    {
        super(WGSMetricsFile.read(metricsFile),
                wgsMetrics -> List.of(
                        new QCValue(QCValueType.REF_COVERAGE_10X, String.valueOf(wgsMetrics.coverage10xPercentage())),
                        new QCValue(QCValueType.REF_COVERAGE_20X, String.valueOf(wgsMetrics.coverage20xPercentage()))));
    }
}
