package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import org.immutables.value.Value;

import net.sf.dynamicreports.report.builder.FieldBuilder;

@Value.Immutable
@Value.Style(allParameters = true)
public abstract class AlterationEvidenceReporterData {
    public static final FieldBuilder<?> SIGNIFICANCE = field("significance", String.class);
    public static final FieldBuilder<String> DRUGS = field("drugs", String.class);
    public static final FieldBuilder<?> SOURCE = field("source", String.class);

    public abstract String getSignificance();

    public abstract String getDrugs();

    public abstract String getSource();

}
