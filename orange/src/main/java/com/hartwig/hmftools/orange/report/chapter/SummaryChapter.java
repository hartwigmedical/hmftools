package com.hartwig.hmftools.orange.report.chapter;

import java.net.MalformedURLException;
import java.text.DecimalFormat;

import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.itextpdf.io.IOException;
import com.itextpdf.io.image.ImageDataFactory;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.property.HorizontalAlignment;

import org.jetbrains.annotations.NotNull;

public class SummaryChapter implements ReportChapter {

    private static final float TABLE_SPACER_HEIGHT = 5;
    private static final DecimalFormat SINGLE_DECIMAL_FORMAT = ReportResources.decimalFormat("#.#");
    private static final DecimalFormat DOUBLE_DECIMAL_FORMAT = ReportResources.decimalFormat("#.##");

    @NotNull
    private final OrangeReport report;

    public SummaryChapter(@NotNull final OrangeReport report) {
        this.report = report;
    }

    @NotNull
    @Override
    public String name() {
        return "ORANGE Report";
    }

    @Override
    public boolean isFullWidth() {
        return false;
    }

    @Override
    public boolean hasCompleteSidebar() {
        return true;
    }

    @Override
    public void render(@NotNull Document reportDocument) {
        reportDocument.add(new Paragraph("Summary").addStyle(ReportResources.chapterTitleStyle()));

        String circosPath = report.plots().purpleCircosPlot();
        try {
            Image circosImage = new Image(ImageDataFactory.create(circosPath));
            circosImage.setMaxHeight(400);
            circosImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
            circosImage.setMarginBottom(8);
            reportDocument.add(circosImage);
        } catch (MalformedURLException e) {
            throw new IOException("Failed to read circos plot image at " + circosPath);
        }
    }
}
