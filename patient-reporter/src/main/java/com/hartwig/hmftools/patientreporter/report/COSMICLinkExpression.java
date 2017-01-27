package com.hartwig.hmftools.patientreporter.report;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.definition.ReportParameters;

class COSMICLinkExpression extends AbstractSimpleExpression<String> {
    public String evaluate(@NotNull final ReportParameters data) {
        return "http://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?id=" + data.getValue(
                PatientDataSource.COSMIC_NR_FIELD.getName());
    }
}
