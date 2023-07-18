package com.hartwig.hmftools.healthchecker.runners;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.metrics.WGSMetrics;
import com.hartwig.hmftools.common.metrics.WGSMetricsFile;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;

public class TumorMetricsChecker extends FileBasedHealthChecker<WGSMetrics>
{
    public TumorMetricsChecker(final String metricsFile) throws IOException
    {
        super(WGSMetricsFile.read(metricsFile),
                wgsMetrics -> List.of(
                        new QCValue(QCValueType.TUM_COVERAGE_30X, String.valueOf(wgsMetrics.coverage30xPercentage())),
                        new QCValue(QCValueType.TUM_COVERAGE_60X, String.valueOf(wgsMetrics.coverage60xPercentage()))));
    }
}
