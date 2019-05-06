package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.common.actionability.ClinicalTrial;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.EvidenceScope;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.Icon;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.cfreport.data.ClinicalTrials;
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

import java.util.List;
import java.util.regex.Pattern;

public class TherapyDetailsChapter implements ReportChapter {

    private final static float COL_WIDTH_DRIVERS = 90;
    private final static float COL_WIDTH_MATCH = 60;
    private final static float COL_WIDTH_LEVEL = 42;
    private final static float COL_WIDTH_RESPONSE_CCMO = 75;
    private final static float COL_WIDTH_SOURCE = 40;
    private final static float COL_WIDTH_TREATMENT_ICONS = 20;
    private final static float COL_WIDTH_TREATMENT_LIST = 180;
    private final static float COL_WIDTH_TREATMENT_LIST_5COL = COL_WIDTH_TREATMENT_LIST + COL_WIDTH_LEVEL;

    private final static String TREATMENT_DELIMITER = " + ";

    private final AnalysedPatientReport patientReport;

    public TherapyDetailsChapter(@NotNull final AnalysedPatientReport patientReport) {
        this.patientReport = patientReport;
    }

    @Override
    public String getName() {
        return "Therapy details";
    }

    @Override
    public final void render(@NotNull final Document reportDocument) {

        Table chapterTable = new Table(1);

        chapterTable.addCell(new Cell()
                .add(createEvidenceTable("Tumor type specific evidence", patientReport.tumorSpecificEvidence()))
                .setPadding(0)
                .setBorder(Border.NO_BORDER));

        chapterTable.addCell(new Cell()
                .add(createClinicalTrialsTable("Clinical trials (NL)", patientReport.clinicalTrials()))
                .setPadding(0)
                .setBorder(Border.NO_BORDER));

        chapterTable.addCell(new Cell()
                .add(createEvidenceTable("Evidence on other tumor types", patientReport.offLabelEvidence()))
                .setPadding(0)
                .setBorder(Border.NO_BORDER));

        // Add legend/footnote that will appear on each page of the chapter
        chapterTable.addFooterCell(new Cell()
                .add(createChapterFootnote())
                .setPadding(0)
                .setBorder(Border.NO_BORDER));

        reportDocument.add(chapterTable);

    }

    @NotNull
    private static Table createEvidenceTable(@NotNull String title, @NotNull final List<EvidenceItem> evidence) {

        // Filter and sort evidence
        final List<EvidenceItem> filteredAndSortedEvidence = EvidenceItems.sort(EvidenceItems.filter(evidence));
        assert(filteredAndSortedEvidence.size() == evidence.size());

        // Handle empty list
        if (filteredAndSortedEvidence.isEmpty()) {
            return TableUtil.createNoneReportTable(title);
        }

        // Create content table
        Table contentTable = TableUtil.createReportContentTable(new float[] {
                COL_WIDTH_DRIVERS,
                COL_WIDTH_MATCH,
                COL_WIDTH_TREATMENT_ICONS,
                COL_WIDTH_TREATMENT_LIST,
                COL_WIDTH_LEVEL,
                COL_WIDTH_RESPONSE_CCMO,
                COL_WIDTH_SOURCE }, new Cell[]  {
                TableUtil.getHeaderCell("Drivers"),
                TableUtil.getHeaderCell("Match"),
                TableUtil.getHeaderCell("Treatments", 2),
                TableUtil.getHeaderCell("Level of evidence"),
                TableUtil.getHeaderCell("Response"),
                TableUtil.getHeaderCell("Source")
        });

        for (EvidenceItem item: filteredAndSortedEvidence) {

            String[] treatments = item.drug().split(Pattern.quote(TREATMENT_DELIMITER));

            contentTable.addCell(TableUtil.getContentCell(item.event()));
            contentTable.addCell(TableUtil.getContentCell(createTreatmentMatchParagraph(item.scope() == EvidenceScope.SPECIFIC)));
            contentTable.addCell(TableUtil.getContentCell(createTreatmentIcons(treatments)).setVerticalAlignment(VerticalAlignment.TOP));
            contentTable.addCell(TableUtil.getContentCell(createTreatmentList(treatments)).setVerticalAlignment(VerticalAlignment.TOP));
            contentTable.addCell(TableUtil.getContentCell(new Paragraph(Icon.createLevelIcon(item.level().readableString()))));
            contentTable.addCell(TableUtil.getContentCell(item.response()));
            contentTable.addCell(TableUtil.getContentCell(new Paragraph(item.source().sourceName()))
                    .setAction(PdfAction.createURI(EvidenceItems.sourceUrl(item))));

        }

        // Create report table that handles page breaks
        return TableUtil.createWrappingReportTable(title, contentTable);

    }

    @NotNull
    private static Table createClinicalTrialsTable(@NotNull String title, @NotNull final List<ClinicalTrial> trials) {

        // Filter and sort trials
        final List<ClinicalTrial> filteredAndSortedTrials = ClinicalTrials.sort(ClinicalTrials.filter(trials));
        assert filteredAndSortedTrials.size() == trials.size();

        // Handle empty list
        if (filteredAndSortedTrials.isEmpty()) {
            return TableUtil.createNoneReportTable(title);
        }

        // Create content table
        final Table contentTable = TableUtil.createReportContentTable(new float[]{COL_WIDTH_DRIVERS,
                        COL_WIDTH_MATCH,
                        COL_WIDTH_TREATMENT_ICONS,
                        COL_WIDTH_TREATMENT_LIST_5COL,
                        COL_WIDTH_RESPONSE_CCMO,
                        COL_WIDTH_SOURCE},
                new Cell[]{
                        TableUtil.getHeaderCell("Drivers"),
                        TableUtil.getHeaderCell("Match"),
                        TableUtil.getHeaderCell("Treatments", 2),
                        TableUtil.getHeaderCell("CCMO"),
                        TableUtil.getHeaderCell("Source")
                });

        for (ClinicalTrial trial: filteredAndSortedTrials) {

            String trialName = trial.acronym();
            contentTable.addCell(TableUtil.getContentCell(trial.event()));
            contentTable.addCell(TableUtil.getContentCell(createTreatmentMatchParagraph(trial.scope() == EvidenceScope.SPECIFIC)));
            contentTable.addCell(TableUtil.getContentCell(createTreatmentIcons(new String[]{trialName})).setVerticalAlignment(VerticalAlignment.TOP));
            contentTable.addCell(TableUtil.getContentCell(trialName).setVerticalAlignment(VerticalAlignment.TOP));
            contentTable.addCell(TableUtil.getContentCell(ClinicalTrials.CCMOId(trial.reference())));
            contentTable.addCell(TableUtil.getContentCell(new Paragraph(trial.source().sourceName())
                    .setAction(PdfAction.createURI(ClinicalTrials.sourceUrl(trial)))));

        }

        // Create report table that handles page breaks
        return TableUtil.createWrappingReportTable(title, contentTable);

    }

    @NotNull
    private static Paragraph createTreatmentMatchParagraph(boolean matchIsSpecific) {
        return new Paragraph()
                .add(Icon.createIcon(matchIsSpecific ? Icon.IconType.MATCH_SPECIFIC : Icon.IconType.MATCH_BROAD))
                .add(new Text(" " + (matchIsSpecific ? "Specific" : "Broad"))
                        .addStyle(ReportResources.tableContentStyle()));
    }

    @NotNull
    private static Paragraph createTreatmentIcons(@NotNull String[] treatments) {
        Paragraph p = new Paragraph();
        for (String treatmentName: treatments) {
            p.add(Icon.createTreatmentIcon(treatmentName.trim()));
        }
        return p;
    }

    @NotNull
    private static Paragraph createTreatmentList(@NotNull String[] treatments) {
        return new Paragraph(String.join(TREATMENT_DELIMITER, treatments))
                .addStyle(ReportResources.tableContentStyle());
    }

    @NotNull
    private static Paragraph createChapterFootnote() {
        return new Paragraph()
                .setKeepTogether(true)
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING)
                .add("The Cancer Genome Interpreter (CGI), OncoKB and CiViC knowledge bases are used to " +
                    "annotate variants of all types with clinical evidence. Only treatment associated evidence with a high " +
                    "level of evidence ( ")
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
