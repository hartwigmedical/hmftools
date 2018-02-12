package com.hartwig.hmftools.patientreporter.report.components;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.definition.ReportParameters;

public class COSMICLinkExpression extends AbstractSimpleExpression<String> {
    @NotNull
    private final FieldBuilder<?> field;

    public COSMICLinkExpression(@NotNull final FieldBuilder<?> field) {
        this.field = field;
    }

    @Override
    public String evaluate(@NotNull final ReportParameters data) {
        return "http://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?id=" + data.getValue(field.getName());
    }
}
