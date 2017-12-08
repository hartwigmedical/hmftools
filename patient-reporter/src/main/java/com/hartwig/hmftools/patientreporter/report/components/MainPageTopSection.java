package com.hartwig.hmftools.patientreporter.report.components;

import static com.hartwig.hmftools.patientreporter.report.Commons.DATE_TIME_FORMAT;
import static com.hartwig.hmftools.patientreporter.report.Commons.dataTableStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.fontStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.tableHeaderStyle;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;

import java.io.IOException;

import com.hartwig.hmftools.patientreporter.SampleReport;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;
import net.sf.dynamicreports.report.constant.VerticalTextAlignment;

public class MainPageTopSection {
    @NotNull
    private static final String REPORT_LOGO_PATH = "pdf/hartwig_logo.jpg";

    @NotNull
    public static ComponentBuilder<?, ?> build(@NotNull final String title, @NotNull final SampleReport report) throws IOException {
        return build(title, report.sampleId(), report.tumorType(), report.tumorPercentageString());
    }

    @NotNull
    public static ComponentBuilder<?, ?> build(@NotNull final String title, @NotNull final String sample, @NotNull final String tumorType,
            @NotNull final String tumorPercentage) {
        // @formatter:off
        final ComponentBuilder<?, ?> mainDiagnosisInfo = cmp.horizontalList(
                cmp.verticalList(
                        cmp.text("Report Date").setStyle(tableHeaderStyle()),
                        cmp.currentDate().setPattern(DATE_TIME_FORMAT).setStyle(dataTableStyle())),
                cmp.verticalList(
                        cmp.text("Primary Tumor Location").setStyle(tableHeaderStyle()),
                        cmp.text(tumorType).setStyle(dataTableStyle())),
                cmp.verticalList(
                        cmp.text("Pathology Tumor Percentage").setStyle(tableHeaderStyle()),
                        cmp.text(tumorPercentage).setStyle(dataTableStyle()))
        );
        return cmp.horizontalList(
                cmp.image(REPORT_LOGO_PATH),
                cmp.verticalList(
                        cmp.text(title + " - " + sample)
                                .setStyle(fontStyle().bold().setFontSize(14)
                                    .setVerticalTextAlignment(VerticalTextAlignment.MIDDLE))
                                .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)
                                .setHeight(50),
                        mainDiagnosisInfo)
        );
        // @formatter:on
    }
}
