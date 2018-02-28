package com.hartwig.hmftools.patientreporter.report.data;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.definition.ReportParameters;

class TranscriptExpression extends AbstractSimpleExpression<String> {

    @NotNull
    private final FieldBuilder<?> transcriptField;

    TranscriptExpression(@NotNull final FieldBuilder<?> transcriptField) {
        this.transcriptField = transcriptField;
    }

    @Override
    public String evaluate(final ReportParameters reportParameters) {
        return "http://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t="
                + reportParameters.getValue(transcriptField.getName());
    }
}
