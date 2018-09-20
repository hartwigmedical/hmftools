package com.hartwig.hmftools.patientreporter.report.components;

import static com.hartwig.hmftools.patientreporter.report.Commons.DATE_TIME_FORMAT;
import static com.hartwig.hmftools.patientreporter.report.Commons.dataTableStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.fontStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.tableHeaderStyle;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;

import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;
import net.sf.dynamicreports.report.constant.VerticalTextAlignment;

public final class MainPageTopSection {

    @NotNull
    private static final String REPORT_LOGO_PATH = "pdf/hartwig_logo.jpg";

    @NotNull
    public static ComponentBuilder<?, ?> buildWithPathologyTumorPercentage(@NotNull final String title,
            @NotNull final SampleReport report) {
        return build(title,
                report.sampleId(),
                report.primaryTumorLocationString(),
                report.cancerSubTypeString(),
                "Pathology Tumor Percentage",
                PatientReportFormat.formatNullablePercent(report.pathologyTumorPercentage()));
    }

    @NotNull
    public static ComponentBuilder<?, ?> buildWithImpliedPurity(@NotNull final String title,
            @NotNull final SampleReport report, @NotNull String purityString) {
        return build(title,
                report.sampleId(),
                report.primaryTumorLocationString(),
                report.cancerSubTypeString(),
                "Implied Tumor Purity",
                purityString);
    }

    @NotNull
    private static ComponentBuilder<?, ?> build(@NotNull String title, @NotNull String sample, @NotNull String primaryTumorLocation,
            @NotNull String cancerSubType, @NotNull String tumorPercentageTitle, @NotNull String tumorPercentage) {
        final ComponentBuilder<?, ?> mainDiagnosisInfo =
                cmp.horizontalList(cmp.verticalList(cmp.text("Report Date").setStyle(tableHeaderStyle()),
                        cmp.currentDate().setPattern(DATE_TIME_FORMAT).setStyle(dataTableStyle())),
                        cmp.verticalList(cmp.text("Primary Tumor Location").setStyle(tableHeaderStyle()),
                                cmp.text(primaryTumorLocation).setStyle(dataTableStyle())),
                        cmp.verticalList(cmp.text("Cancer Subtype").setStyle(tableHeaderStyle()),
                                cmp.text(cancerSubType).setStyle(dataTableStyle())),
                        cmp.verticalList(cmp.text(tumorPercentageTitle).setStyle(tableHeaderStyle()),
                                cmp.text(tumorPercentage).setStyle(dataTableStyle()))).setHeight(32);
        return cmp.verticalList(cmp.horizontalList(cmp.image(REPORT_LOGO_PATH).setWidth(42).setHeight(65),
                cmp.text(title + " - " + sample)
                        .setStyle(fontStyle().bold().setFontSize(14).setVerticalTextAlignment(VerticalTextAlignment.MIDDLE))
                        .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER),
                cmp.horizontalGap(40)), cmp.verticalGap(3), mainDiagnosisInfo);
    }
}
