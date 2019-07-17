package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import java.util.List;
import java.util.regex.Pattern;

import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.EvidenceItemMerger;
import com.hartwig.hmftools.common.actionability.EvidenceScope;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.Icon;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.cfreport.data.EvidenceDrugTypeMerger;
import com.hartwig.hmftools.patientreporter.cfreport.data.EvidenceItems;
import com.itextpdf.kernel.pdf.action.PdfAction;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.borders.Border;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.element.Text;
import com.itextpdf.layout.property.VerticalAlignment;

import org.jetbrains.annotations.NotNull;

public class TherapyDetailsChapterOffLabel implements ReportChapter {

    private final static float COL_WIDTH_DRIVERS = 110;
    private final static float COL_WIDTH_MATCH = 60;
    private final static float COL_WIDTH_LEVEL = 42;
    private final static float COL_WIDTH_RESPONSE_CCMO = 75;
    private final static float COL_WIDTH_SOURCE = 40;
    private final static float COL_WIDTH_TREATMENT_ICONS = 25;
    private final static float COL_WIDTH_TREATMENT_LIST = 180;

    private final static String TREATMENT_DELIMITER = " + ";

    @NotNull
    private final AnalysedPatientReport patientReport;

    public TherapyDetailsChapterOffLabel(@NotNull final AnalysedPatientReport patientReport) {
        this.patientReport = patientReport;
    }

    @NotNull
    @Override
    public String name() {
        return "Therapy details (Other tumor types)";
    }

    @Override
    public final void render(@NotNull final Document reportDocument) {
        Table chapterTable = new Table(1);

        chapterTable.addCell(new Cell().add(createEvidenceTable("Evidence on other tumor types", patientReport.offLabelEvidence()))
                .setPadding(0)
                .setBorder(Border.NO_BORDER));

        chapterTable.addFooterCell(new Cell().add(createChapterFootnote()).setPadding(0).setBorder(Border.NO_BORDER));

        reportDocument.add(chapterTable);
    }

    @NotNull
    private static Table createEvidenceTable(@NotNull String title, @NotNull List<EvidenceItem> evidence) {
        if (evidence.isEmpty()) {
            return TableUtil.createNoneReportTable(title);
        }

        Table contentTable = TableUtil.createReportContentTable(new float[] { COL_WIDTH_DRIVERS, COL_WIDTH_MATCH, COL_WIDTH_TREATMENT_ICONS,
                        COL_WIDTH_TREATMENT_LIST, COL_WIDTH_LEVEL, COL_WIDTH_RESPONSE_CCMO, COL_WIDTH_SOURCE },
                new Cell[] { TableUtil.createHeaderCell("Variant"), TableUtil.createHeaderCell("Match"),
                        TableUtil.createHeaderCell("Treatment", 2), TableUtil.createHeaderCell("Level of evidence"),
                        TableUtil.createHeaderCell("Response"), TableUtil.createHeaderCell("Source") });

      //  List<EvidenceItemMerger> mergedItems = EvidenceDrugTypeMerger.merge(evidence);

        final List<EvidenceItem> sortedEvidence = EvidenceItems.sort(evidence);
        for (EvidenceItem item : sortedEvidence) {
            String[] treatments = item.drug().split(Pattern.quote(TREATMENT_DELIMITER));

            contentTable.addCell(TableUtil.createContentCell(item.event()));
            contentTable.addCell(TableUtil.createContentCell(createTreatmentMatchParagraph(item.scope() == EvidenceScope.SPECIFIC)));
            contentTable.addCell(TableUtil.createContentCell(createTreatmentIcons(treatments)).setVerticalAlignment(VerticalAlignment.TOP));
            contentTable.addCell(TableUtil.createContentCell(createTreatmentList(treatments)).setVerticalAlignment(VerticalAlignment.TOP));
            contentTable.addCell(TableUtil.createContentCell(new Paragraph(Icon.createLevelIcon(item.level().readableString()))));
            contentTable.addCell(TableUtil.createContentCell(item.response()));
            contentTable.addCell(TableUtil.createContentCell(new Paragraph(item.source()
                    .sourceName()).addStyle(ReportResources.dataHighlightLinksStyle()))
                    .setAction(PdfAction.createURI(EvidenceItems.sourceUrl(item))));
        }

        return TableUtil.createWrappingReportTable(title, contentTable);
    }

    @NotNull
    private static Paragraph createTreatmentMatchParagraph(boolean matchIsSpecific) {
        return new Paragraph().add(Icon.createIcon(matchIsSpecific ? Icon.IconType.MATCH_SPECIFIC : Icon.IconType.MATCH_BROAD))
                .add(new Text(" " + (matchIsSpecific ? "Specific" : "Broad")).addStyle(ReportResources.tableContentStyle()));
    }

    @NotNull
    private static Paragraph createTreatmentIcons(@NotNull String[] treatments) {
        Paragraph p = new Paragraph();
        for (String treatmentName : treatments) {
            p.add(Icon.createTreatmentIcon(treatmentName.trim()));
        }
        return p;
    }

    @NotNull
    private static Paragraph createTreatmentList(@NotNull String[] treatments) {
        return new Paragraph(String.join(TREATMENT_DELIMITER, treatments)).addStyle(ReportResources.tableContentStyle());
    }

    @NotNull
    private static Paragraph createChapterFootnote() {
        return new Paragraph().setKeepTogether(true)
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING)
                .add("The Cancer Genome Interpreter (CGI), OncoKB and CiViC knowledge bases are used to "
                        + "annotate variants of all types with clinical evidence. Only treatment associated evidence with a high "
                        + "level of evidence ( ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_A))
                .add(" validated association; ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_B))
                .add(" strong clinical evidence) are reported here. Potential evidence items with a lower level of evidence ( ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_C))
                .add(" case study, limited clinical evidence; ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_D))
                .add(" pre-clinical) are not reported.")
                .addStyle(ReportResources.subTextStyle());
    }
}
