package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.definition.ReportParameters;

@Value.Immutable
@Value.Style(allParameters = true)
public abstract class AlterationEvidence {
    public static final FieldBuilder<?> SIGNIFICANCE = field("significance", String.class);
    public static final FieldBuilder<String> DRUGS = field("drugs", String.class);
    public static final FieldBuilder<?> SOURCE = field("source", String.class);
    public static final FieldBuilder<?> SOURCE_URL = field("sourceUrl", String.class);

    public abstract String getSignificance();

    public abstract String getDrugs();

    public abstract String getSource();

    public abstract String getSourceUrl();

    @NotNull
    public static AbstractSimpleExpression<String> sourceHyperlink() {
        return new AbstractSimpleExpression<String>() {
            @Override
            public String evaluate(@NotNull final ReportParameters data) {
                return data.getValue(SOURCE_URL.getName());
            }
        };
    }
}
