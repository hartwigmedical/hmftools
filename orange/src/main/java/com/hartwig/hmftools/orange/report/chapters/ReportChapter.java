package com.hartwig.hmftools.orange.report.chapters;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;

public interface ReportChapter
{
    String name();

    boolean isLandscape();

    JasperReportBuilder buildReport();
}
