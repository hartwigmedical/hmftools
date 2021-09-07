package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import java.util.List;

import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.patientreporter.algo.GenomicAnalysis;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.cfreport.data.ClinicalTrials;
import com.itextpdf.kernel.pdf.action.PdfAction;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.borders.Border;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.VerticalAlignment;

import org.jetbrains.annotations.NotNull;

public class TherapyDetailsChapterOnLabel implements ReportChapter {

    private static final float COL_WIDTH_GENOMIC_EVENT = 110;
    private static final float COL_WIDTH_TREATMENT_ICONS = 25;
    private static final float COL_WIDTH_TRIAL_NAME = 222;
    private static final float COL_WIDTH_SOURCE = 40;

    @NotNull
    private final GenomicAnalysis genomicAnalysis;

    public TherapyDetailsChapterOnLabel(@NotNull final GenomicAnalysis genomicAnalysis) {
        this.genomicAnalysis = genomicAnalysis;
    }

    @NotNull
    @Override
    public String name() {
        return "Therapy details (Tumor type specific)";
    }

    @Override
    public void render(@NotNull Document reportDocument) {
        Table chapterTable = new Table(1);

        chapterTable.addCell(new Cell().add(TherapyDetailsChapterFunctions.createEvidenceTable("Tumor type specific evidence",
                genomicAnalysis.tumorSpecificEvidence())).setPadding(0).setBorder(Border.NO_BORDER));

        chapterTable.addCell(new Cell().add(createClinicalTrialsTable(genomicAnalysis.clinicalTrials()))
                .setPadding(0)
                .setBorder(Border.NO_BORDER));

        chapterTable.addFooterCell(new Cell().add(TherapyDetailsChapterFunctions.createChapterFootnoteOnLabel())
                .setPadding(0)
                .setBorder(Border.NO_BORDER));

        reportDocument.add(chapterTable);
    }

    @NotNull
    private static Table createClinicalTrialsTable(@NotNull List<ProtectEvidence> trials) {
        String title = "Tumor type specific clinical trials (NL)";

        if (trials.isEmpty()) {
            return TableUtil.createNoneReportTable(title);
        }

        Table contentTable =
                TableUtil.createReportContentTable(new float[] { COL_WIDTH_GENOMIC_EVENT, COL_WIDTH_TREATMENT_ICONS, COL_WIDTH_TRIAL_NAME,
                                COL_WIDTH_SOURCE },
                        new Cell[] { TableUtil.createHeaderCell("Genomic event"), TableUtil.createHeaderCell("Trial", 2),
                                TableUtil.createHeaderCell("Source") });

        for (ProtectEvidence trial : ClinicalTrials.sort(trials)) {
            String trialName = trial.treatment();
            contentTable.addCell(TableUtil.createContentCell(trial.genomicEvent()));
            contentTable.addCell(TableUtil.createContentCell(TherapyDetailsChapterFunctions.createTreatmentIcons(trialName))
                    .setVerticalAlignment(VerticalAlignment.TOP));
            contentTable.addCell(TableUtil.createContentCell(trialName).setVerticalAlignment(VerticalAlignment.TOP));
            contentTable.addCell(TableUtil.createContentCell(new Paragraph(ClinicalTrials.source(trial)).addStyle(ReportResources.dataHighlightLinksStyle()))
                    .setAction(PdfAction.createURI(ClinicalTrials.trialUrl(trial))));
        }

        contentTable.addCell(TableUtil.createLayoutCell(1, contentTable.getNumberOfColumns())
                .setPaddingTop(10)
                .add(new Paragraph("Potential eligibility for DRUP is dependent on tumor type details therefore certain tumor types "
                        + "may not be eligible for the DRUP.").addStyle(ReportResources.subTextStyle())));
        return TableUtil.createWrappingReportTable(title, contentTable);
    }
}
