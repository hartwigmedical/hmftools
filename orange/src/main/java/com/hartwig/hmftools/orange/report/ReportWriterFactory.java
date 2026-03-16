package com.hartwig.hmftools.orange.report;

import com.hartwig.hmftools.orange.OrangeConfig;

public final class ReportWriterFactory
{
    public static ReportWriter createToDiskWriter(final OrangeConfig config)
    {
        return new ReportWriter(true, config, new PlotPathResolver(config.OutputDir));
    }

    public static ReportWriter createInMemoryWriter()
    {
        return new ReportWriter(false, null, new PlotPathResolver(null));
    }
}
