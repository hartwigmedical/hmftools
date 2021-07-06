package com.hartwig.hmftools.orange.report.chapters;

import java.net.MalformedURLException;

import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.components.TableUtil;
import com.itextpdf.io.IOException;
import com.itextpdf.io.image.ImageDataFactory;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.HorizontalAlignment;

import org.jetbrains.annotations.NotNull;

public class FrontPageChapter implements ReportChapter {

    @NotNull
    private final OrangeReport report;

    public FrontPageChapter(@NotNull final OrangeReport report) {
        this.report = report;
    }

    @NotNull
    @Override
    public String name() {
        return "Front Page";
    }

    @Override
    public void render(@NotNull Document document) {
        Table content = TableUtil.createReportContentTable(new float[] { 1, 1 },
                new Cell[] { TableUtil.createHeaderCell("Clinical Tumor Type"), TableUtil.createHeaderCell("Molecular Tumor Type") });

        content.addCell(TableUtil.createContentCell("8923 - Melanoma"));
        content.addCell(TableUtil.createContentCell("Melanoma"));

        document.add(TableUtil.createWrappingReportTable(content));

        String circosPath = report.plots().purpleCircosPlot();
        try {
            Image circosImage = new Image(ImageDataFactory.create(circosPath));
            circosImage.setMaxHeight(400);
            circosImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
            circosImage.setMarginBottom(8);
            document.add(circosImage);
        } catch (MalformedURLException e) {
            throw new IOException("Failed to read circos plot image at " + circosPath);
        }
    }
}
