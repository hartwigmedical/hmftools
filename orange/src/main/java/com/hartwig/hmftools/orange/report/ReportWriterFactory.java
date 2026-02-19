package com.hartwig.hmftools.orange.report;

import com.hartwig.hmftools.orange.OrangeConfig;

public final class ReportWriterFactory
{
    public static ReportWriter createToDiskWriter(final OrangeConfig config)
    {
        return new ReportWriter(true, config.OutputDir, new PlotPathResolver(config.OutputDir), config.AddDisclaimer);
    }

    public static ReportWriter createInMemoryWriter()
    {
        return new ReportWriter(false, null, new PlotPathResolver(null), true);
    }
}
