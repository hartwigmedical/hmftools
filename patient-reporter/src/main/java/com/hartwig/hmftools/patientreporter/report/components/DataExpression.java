package com.hartwig.hmftools.patientreporter.report.components;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.definition.ReportParameters;

public class DataExpression extends AbstractSimpleExpression<String> {

    @NotNull
    private final String field;

    @NotNull
    public static DataExpression fromField(@NotNull FieldBuilder<?> builder) {
        return new DataExpression(builder.getName());
    }

    private DataExpression(@NotNull final String field) {
        this.field = field;
    }

    @Override
    @NotNull
    public String evaluate(@NotNull final ReportParameters reportParameters) {
        return reportParameters.getValue(field);
    }
}
