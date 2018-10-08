package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

public class ActionabilityVariantsDataSource {
    public static final FieldBuilder<?> EVENT = field("event", String.class);
    public static final FieldBuilder<?> MATCHING_CANCERTYPE = field("matching cancerType", String.class);
    public static final FieldBuilder<?> SOURCE = field("source", String.class);
    public static final FieldBuilder<?> DRUG = field("drug", String.class);
    public static final FieldBuilder<?> DRUGS_TYPE = field("drugs type", String.class);
    public static final FieldBuilder<?> LEVEL = field("level", String.class);
    public static final FieldBuilder<?> RESPONSE = field("response", String.class);


    private ActionabilityVariantsDataSource() {
    }

    @NotNull
    public static JRDataSource fromActionabilityVariants() {
        final DRDataSource actionabilityVariantsDatasource = new DRDataSource(EVENT.getName(),
                MATCHING_CANCERTYPE.getName(),
                SOURCE.getName(),
                DRUG.getName(),
                DRUGS_TYPE.getName(),
                LEVEL.getName(),
                RESPONSE.getName());




        return actionabilityVariantsDatasource;
    }

    @NotNull
    public static FieldBuilder<?>[] actionabilityFields() {
        return new FieldBuilder<?>[] { EVENT, MATCHING_CANCERTYPE, SOURCE, DRUG, DRUGS_TYPE, LEVEL, RESPONSE};
    }
}
