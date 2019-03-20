package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.BodyText;
import com.hartwig.hmftools.patientreporter.cfreport.components.LineDivider;
import com.hartwig.hmftools.patientreporter.cfreport.components.SectionTitle;
import com.hartwig.hmftools.patientreporter.cfreport.components.tables.TableHelper;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.Style;
import com.itextpdf.layout.borders.Border;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;
import com.itextpdf.layout.property.VerticalAlignment;
import com.sun.xml.internal.ws.util.StringUtils;
import com.sun.xml.internal.ws.wsdl.writer.document.soap.Body;
import org.jcp.xml.dsig.internal.dom.ApacheCanonicalizer;
import org.jetbrains.annotations.NotNull;

import java.util.StringJoiner;

public class SummaryChapter extends ReportChapter {

    private final static float TABLE_SPACER_HEIGHT = 10;

    @Override
    public final String getName() {
        return "Summary";
    }

    @Override
    public final ChapterType getChapterType() {
        return ChapterType.SummaryChapter;
    }

    @Override
    protected final void renderChapterContent(@NotNull Document report) {
        renderTumorLocationAndType(report);
        renderSummaryText(report);
        renderTreatmentIndications(report);
        renderTumorCharacteristics(report);
        renderGenomicAlterations(report);
    }

    private final void renderTumorLocationAndType(@NotNull Document report) {

        Div div = new Div();
        div.setKeepTogether(true);
        div.setWidth(getContentWidth());

    }

    private final void renderSummaryText(@NotNull Document report) {

        Div div = initializeSummarySectionDiv(getContentWidth());
        div.add(SectionTitle.getSectionTitle("Summary"));

        // Add content to div
        String content = "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Vivamus eget porta turpis. Lorem ipsum dolor sit amet, " +
                "consectetur adipiscing elit. Nullam interdum sodales ullamcorper. Nulla vestibulum ipsum quis ipsum congue, quis commodo velit " +
                "condimentum. Suspendisse eget nulla egestas, fermentum urna ut, bibendum ipsum. Nulla varius, dui elementum faucibus ultricies, " +
                "nisi velit dignissim arcu, nec feugiat dui magna eu felis. Maecenas at odio pharetra, sodales velit vitae, gravida mauris. Pellentesque " +
                "id ultrices diam. Integer non ex ut neque auctor pellentesque. Ut et nibh faucibus, pretium erat efficitur, vehicula lorem.";
        Paragraph p = BodyText.getParagraph(content);
        p.setWidth(getContentWidth());
        div.add(p);

        report.add(div);

    }

    private final void renderTreatmentIndications(@NotNull Document report) {

        // Initialize div
        Div div = initializeSummarySectionDiv(getContentWidth());

        // Initialize table
        Table table = new Table(UnitValue.createPercentArray(new float[] {1, 1}));
        table.setWidth(getContentWidth());

        table.addCell(TableHelper.getLayoutCell().add(SectionTitle.getSectionTitle("Treatment indications")));
        table.addCell(TableHelper.getLayoutCell().add(BodyText.getParagraph("Summary of number of alterations with number of treatment indication and/or clinical studies")));
        table.addCell(TableHelper.getLayoutCell(1, 2).setHeight(TABLE_SPACER_HEIGHT)); // Spacer

        table.addCell(getBottomAlignedLayoutCell()
                .add(BodyText.getParagraph("Gene alteration(s) with therapy indication(s)")));
        table.addCell(getBottomAlignedLayoutCell()
                .add(new Paragraph(String.format("%d (%d treatments)", 1, 8))
                .addStyle(ReportResources.dataHighlightStyle())));

        table.addCell(getBottomAlignedLayoutCell()
                .add(BodyText.getParagraph("Gene alteration(s) with clinical study eligibility")));
        table.addCell(getBottomAlignedLayoutCell()
                .add(new Paragraph(String.format("%d (%d studies)", 1, 8))
                .addStyle(ReportResources.dataHighlightStyle())));

        div.add(table);

        report.add(div);

    }

    private final void renderTumorCharacteristics(@NotNull Document report) {

        // Initialize div
        Div div = initializeSummarySectionDiv(getContentWidth());

        // Initialize table
        Table table = new Table(UnitValue.createPercentArray(new float[] {1, .33f, .66f}));
        table.setWidth(getContentWidth());

        table.addCell(TableHelper.getLayoutCell().add(SectionTitle.getSectionTitle("Tumor characteristics summary")));
        table.addCell(TableHelper.getLayoutCell(1, 2).add(BodyText.getParagraph("Whole genome sequencing based tumor characteristics.")));
        table.addCell(TableHelper.getLayoutCell(1, 3).setHeight(TABLE_SPACER_HEIGHT)); // Spacer

        table.addCell(getBottomAlignedLayoutCell().add(BodyText.getParagraph("Tumor purity of biopsy")));
        table.addCell(getBottomAlignedLayoutCell().add(new Paragraph("[VAL]")));
        table.addCell(getBottomAlignedLayoutCell().add(new Paragraph("[BAR]")));

        table.addCell(getBottomAlignedLayoutCell().add(BodyText.getParagraph("Average tumor ploidy")));
        table.addCell(getBottomAlignedLayoutCell().add(new Paragraph("[VAL]")));
        table.addCell(getBottomAlignedLayoutCell().add(new Paragraph("[BAR]")));

        table.addCell(getBottomAlignedLayoutCell().add(BodyText.getParagraph("Tumor mutational load")));
        table.addCell(getBottomAlignedLayoutCell().add(new Paragraph("[VAL]")));
        table.addCell(getBottomAlignedLayoutCell().add(new Paragraph("[BAR]")));

        table.addCell(getBottomAlignedLayoutCell().add(BodyText.getParagraph("Microsatellite (in)stability")));
        table.addCell(getBottomAlignedLayoutCell().add(new Paragraph("[VAL]")));
        table.addCell(getBottomAlignedLayoutCell().add(new Paragraph("[BAR]")));

        div.add(table);

        report.add(div);

    }

    private final void renderGenomicAlterations(@NotNull Document report) {

        // Initialize div
        Div div = initializeSummarySectionDiv(getContentWidth());

        // Initialize table
        Table table = new Table(UnitValue.createPercentArray(new float[] {1, 1}));
        table.setWidth(getContentWidth());
        table.addCell(TableHelper.getLayoutCell().add(SectionTitle.getSectionTitle("Genomic alterations\n summary")));
        table.addCell(TableHelper.getLayoutCell().add(BodyText.getParagraph("Summary on genomic alterations " +
                "(somatic variants, copy number changes, gene disruptions and gene fusions).")));
        table.addCell(TableHelper.getLayoutCell(1, 2).setHeight(TABLE_SPACER_HEIGHT)); // Spacer

        // Genes with driver variant
        String[] driverVariantGenes = {"CDKN2A", "BRAF"};
        table.addCell(getBottomAlignedLayoutCell().add(BodyText.getParagraph("Genes with driver variant")));
        table.addCell(getGeneListCell(driverVariantGenes));

        // Reported variants
        int reportedVariants = 4;
        Style reportedVariantsStyle = (reportedVariants > 0) ? ReportResources.dataHighlightStyle() : ReportResources.dataHighlightNaStyle();
        table.addCell(getBottomAlignedLayoutCell().add(BodyText.getParagraph("Nr. of reported variants")));
        table.addCell(getBottomAlignedLayoutCell().add(new Paragraph(String.valueOf(reportedVariants))
                .addStyle(reportedVariantsStyle)));

        // Copy gain genes
        String[] copyGainGenes = {};
        table.addCell(getBottomAlignedLayoutCell().add(BodyText.getParagraph("Genes with copy-gain")));
        table.addCell(getGeneListCell(copyGainGenes));

        // Copy loss genes
        String[] copyLossGenes = {"PTEN"};
        table.addCell(getBottomAlignedLayoutCell().add(BodyText.getParagraph("Genes with copy-loss")));
        table.addCell(getGeneListCell(copyLossGenes));

        // Gene fusions
        String[] fusionGenes = {};
        table.addCell(getBottomAlignedLayoutCell().add(BodyText.getParagraph("Gene fusions")));
        table.addCell(getGeneListCell(fusionGenes));

        div.add(table);

        report.add(div);


    }

    @NotNull
    private static final Div initializeSummarySectionDiv(float width) {

        Div div = new Div();
        div.setKeepTogether(true);
        div.setWidth(width);

        // Add divider and section title
        div.add(LineDivider.getLineDivider(width));

        return div;

    }

    @NotNull
    private static final Cell getBottomAlignedLayoutCell() {
        Cell c = TableHelper.getLayoutCell()
                .setVerticalAlignment(VerticalAlignment.BOTTOM);
        return c;
    }

    @NotNull
    private static final Cell getGeneListCell(@NotNull String[] genes) {

        // Concatenate genes
        String geneString;
        if (genes.length == 0) {
            geneString = "NONE";
        } else {

            StringJoiner joiner = new StringJoiner(", ");
            for (String s: genes) {
                joiner.add(s);
            }
            geneString = joiner.toString();

        }

        // Fetch style
        Style style = genes.length > 0 ? ReportResources.dataHighlightStyle() : ReportResources.dataHighlightNaStyle();

        // Build table
        Cell c = getBottomAlignedLayoutCell()
                .add(new Paragraph(geneString))
                .addStyle(style);
                return c;

    }


}
