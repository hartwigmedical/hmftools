package com.hartwig.hmftools.patientreporter.report.components;

import static com.hartwig.hmftools.patientreporter.report.Commons.DATE_TIME_FORMAT;
import static com.hartwig.hmftools.patientreporter.report.Commons.dataTableStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.fontStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.tableHeaderStyle;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;

import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.SampleReport;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;
import net.sf.dynamicreports.report.constant.VerticalTextAlignment;

public final class MainPageTopSection {

    @NotNull
    public static ComponentBuilder<?, ?> build(@NotNull String title, @NotNull PatientReport report) {
        SampleReport sampleReport = report.sampleReport();
        return build(sampleReport.buildReportTitle(title),
                sampleReport.primaryTumorLocationString(),
                sampleReport.cancerSubTypeString(),
                report.logoCompanyPath());
    }

    @NotNull
    private static ComponentBuilder<?, ?> build(@NotNull String title, @NotNull String primaryTumorLocation, @NotNull String cancerSubType,
            @NotNull String companyLogoPath) {
        final ComponentBuilder<?, ?> mainDiagnosisInfo =
                cmp.horizontalList(cmp.verticalList(cmp.text("Report Date").setStyle(tableHeaderStyle().setPadding(2)),
                        cmp.currentDate().setPattern(DATE_TIME_FORMAT).setStyle(dataTableStyle().setPadding(2))),
                        cmp.verticalList(cmp.text("Primary Tumor Location").setStyle(tableHeaderStyle().setPadding(2)),
                                cmp.text(primaryTumorLocation).setStyle(dataTableStyle().setPadding(2))),
                        cmp.verticalList(cmp.text("Cancer Subtype").setStyle(tableHeaderStyle().setPadding(2)),
                                cmp.text(cancerSubType).setStyle(dataTableStyle().setPadding(2))));

        return cmp.verticalList(cmp.horizontalList(cmp.image(companyLogoPath).setWidth(42).setHeight(65),
                cmp.text(title)
                        .setStyle(fontStyle().bold().setFontSize(14).setVerticalTextAlignment(VerticalTextAlignment.MIDDLE))
                        .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER),
                cmp.horizontalGap(40)), cmp.verticalGap(3), mainDiagnosisInfo);
    }
}
