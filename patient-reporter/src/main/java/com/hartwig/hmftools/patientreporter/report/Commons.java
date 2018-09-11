package com.hartwig.hmftools.patientreporter.report;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.stl;

import java.awt.Color;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;

import com.hartwig.hmftools.patientreporter.PatientReporterApplication;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.component.VerticalListBuilder;
import net.sf.dynamicreports.report.builder.style.StyleBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;
import net.sf.dynamicreports.report.constant.VerticalTextAlignment;

public final class Commons {

    public static final String TITLE_SEQUENCE = "HMF Sequencing Report v" + PatientReporterApplication.VERSION;

    public static final String DATE_TIME_FORMAT = "dd-MMM-yyyy";
    public static final String HARTWIG_ADDRESS = "Hartwig Medical Foundation, Science Park 408, 1098XH Amsterdam";

    public static final int SECTION_VERTICAL_GAP = 25;
    public static final int HEADER_TO_DETAIL_VERTICAL_GAP = 8;

    private static final String FONT = "Tinos";
    private static final String MONOSPACE_FONT = "Inconsolata";
    private static final Color BORKIE_COLOR = new Color(221, 235, 247);
    private static final int TABLE_PADDING = 1;

    private static final int TEXT_HEADER_INDENT = 30;
    private static final int TEXT_DETAIL_INDENT = 40;
    private static final int TEXT_END_OF_LINE_GAP = 30;
    private static final int LIST_INDENT = 5;
    private static final int DETAIL_TO_DETAIL_VERTICAL_GAP = 4;

    private Commons() {
    }

    @NotNull
    public static VerticalListBuilder toList(@NotNull final String title, @NotNull final Iterable<String> lines) {
        final VerticalListBuilder list = cmp.verticalList();
        list.add(cmp.horizontalList(cmp.horizontalGap(TEXT_HEADER_INDENT), cmp.text(title).setStyle(fontStyle().bold().setFontSize(11))),
                cmp.verticalGap(HEADER_TO_DETAIL_VERTICAL_GAP));
        boolean isFirst = true;
        for (final String line : lines) {
            if (!isFirst) {
                list.add(cmp.verticalGap(DETAIL_TO_DETAIL_VERTICAL_GAP));
            }
            list.add(cmp.horizontalList(cmp.horizontalGap(TEXT_DETAIL_INDENT),
                    cmp.text("- ").setStyle(fontStyle()).setWidth(LIST_INDENT),
                    cmp.text(line).setStyle(fontStyle()),
                    cmp.horizontalGap(TEXT_END_OF_LINE_GAP)));

            isFirst = false;
        }
        return list;
    }

    @NotNull
    public static StyleBuilder tableHeaderStyle() {
        return fontStyle().bold()
                .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)
                .setVerticalTextAlignment(VerticalTextAlignment.MIDDLE)
                .setFontSize(10)
                .setBorder(stl.pen1Point())
                .setBackgroundColor(BORKIE_COLOR)
                .setPadding(TABLE_PADDING);
    }

    @NotNull
    public static StyleBuilder dataStyle() {
        return fontStyle().setFontSize(8)
                .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)
                .setVerticalTextAlignment(VerticalTextAlignment.MIDDLE)
                .setPadding(TABLE_PADDING);
    }

    @NotNull
    public static StyleBuilder dataTableStyle() {
        return dataStyle().setBorder(stl.pen1Point());
    }

    @NotNull
    public static StyleBuilder fontStyle() {
        return stl.style().setFontName(FONT);
    }

    @NotNull
    public static StyleBuilder monospaceFontStyle() {
        return stl.style().setFontName(MONOSPACE_FONT);
    }

    @NotNull
    public static JasperReportBuilder baseTable() {
        return report().setColumnStyle(dataStyle()).setColumnTitleStyle(tableHeaderStyle()).highlightDetailEvenRows();
    }

    @NotNull
    public static JasperReportBuilder monospaceBaseTable() {
        return report().setColumnStyle(dataStyle().setFontName(MONOSPACE_FONT))
                .setColumnTitleStyle(tableHeaderStyle())
                .highlightDetailEvenRows();
    }

    @NotNull
    public static StyleBuilder sectionHeaderStyle() {
        return fontStyle().bold().setFontSize(12).setHorizontalTextAlignment(HorizontalTextAlignment.CENTER);
    }

    @NotNull
    public static StyleBuilder linkStyle() {
        return dataStyle().setForegroundColor(Color.BLUE);
    }

    @NotNull
    public static String formattedDate(@Nullable final LocalDate date) {
        final DateTimeFormatter formatter = DateTimeFormatter.ofPattern(DATE_TIME_FORMAT);
        return date != null ? formatter.format(date) : "?";
    }
}
