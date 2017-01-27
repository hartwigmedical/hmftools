package com.hartwig.hmftools.patientreporter.report;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.definition.ReportParameters;

class TranscriptLinkExpression extends AbstractSimpleExpression<String> {
    public String evaluate(@NotNull final ReportParameters data) {
        return "http://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=" + data.getValue(
                PatientDataSource.TRANSCRIPT_FIELD.getName());
    }
}
