package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;

import org.immutables.value.Value;

import net.sf.dynamicreports.report.builder.FieldBuilder;

@Value.Immutable
@Value.Style(allParameters = true)
public abstract class AlterationReporterData {
    public static final FieldBuilder<?> VARIANT = field("variant", String.class);
    public static final FieldBuilder<?> PREDICTED_EFFECT = field("predictedEffect", String.class);
    public static final FieldBuilder<?> ALTERATION = field("alteration", String.class);
    public static final FieldBuilder<?> LOH = field("loh", String.class);
    public static final FieldBuilder<?> SUBCLONAL = field("subclonal", String.class);

    public abstract String getVariant();

    public abstract String getPredictedEffect();

    public abstract String getLoh();

    public abstract String getSubclonal();

    public abstract List<AlterationEvidenceReporterData> getEvidence();

    public String getAlteration() {
        return getVariant() + "\n" + getPredictedEffect();
    }
}
