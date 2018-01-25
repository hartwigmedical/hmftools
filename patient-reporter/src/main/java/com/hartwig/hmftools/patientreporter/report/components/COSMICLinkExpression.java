package com.hartwig.hmftools.patientreporter.report.components;

import com.hartwig.hmftools.patientreporter.report.data.VariantDataSource;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.definition.ReportParameters;

public class COSMICLinkExpression extends AbstractSimpleExpression<String> {
    @Override
    public String evaluate(@NotNull final ReportParameters data) {
        return "http://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?id=" + data.getValue(
                VariantDataSource.COSMIC_NR_FIELD.getName());
    }
}
