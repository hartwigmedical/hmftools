package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import java.util.List;
import java.util.regex.Pattern;

import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.Icon;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.cfreport.data.EvidenceItems;
import com.itextpdf.kernel.pdf.action.PdfAction;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.VerticalAlignment;

import org.jetbrains.annotations.NotNull;

final class TherapyDetailsChapterFunctions {

    private static final float COL_WIDTH_GENOMIC_EVENT = 110;
    private static final float COL_WIDTH_DRUG_ICONS = 25;
    private static final float COL_WIDTH_DRUG_LIST = 180;
    private static final float COL_WIDTH_LEVEL = 42;
    private static final float COL_WIDTH_RESPONSE = 75;
    private static final float COL_WIDTH_SOURCE = 40;

    private static final String TREATMENT_DELIMITER = " + ";

    private TherapyDetailsChapterFunctions() {
    }

    @NotNull
    static Table createEvidenceTable(@NotNull String title, @NotNull List<ProtectEvidence> evidence) {
        if (evidence.isEmpty()) {
            return TableUtil.createNoneReportTable(title);
        }

        Table contentTable =
                TableUtil.createReportContentTable(new float[] { COL_WIDTH_GENOMIC_EVENT, COL_WIDTH_DRUG_ICONS, COL_WIDTH_DRUG_LIST,
                                COL_WIDTH_LEVEL, COL_WIDTH_RESPONSE, COL_WIDTH_SOURCE },
                        new Cell[] { TableUtil.createHeaderCell("Genomic event"), TableUtil.createHeaderCell("Treatment", 2),
                                TableUtil.createHeaderCell("Level"), TableUtil.createHeaderCell("Direction"),
                                TableUtil.createHeaderCell("Sources") });

        for (ProtectEvidence item : EvidenceItems.sort(evidence)) {
            contentTable.addCell(TableUtil.createContentCell(item.genomicEvent()));
            contentTable.addCell(TableUtil.createContentCell(createTreatmentIcons(item.treatment()).setVerticalAlignment(VerticalAlignment.TOP)));
            contentTable.addCell(TableUtil.createContentCell(new Paragraph(item.treatment()).addStyle(ReportResources.tableContentStyle())
                    .setVerticalAlignment(VerticalAlignment.TOP)));
            contentTable.addCell(TableUtil.createContentCell(new Paragraph(Icon.createLevelIcon(item.level().name()))));
            contentTable.addCell(TableUtil.createContentCell(item.direction().name()));
            contentTable.addCell(TableUtil.createContentCell(EvidenceItems.evidenceUrl(item).isEmpty()
                    ? new Paragraph(EvidenceItems.sources(item))
                    : new Paragraph(EvidenceItems.sources(item)).addStyle(ReportResources.dataHighlightLinksStyle()))
                    .setAction(PdfAction.createURI(EvidenceItems.evidenceUrl(item))));
        }

        return TableUtil.createWrappingReportTable(title, contentTable);
    }

    @NotNull
    static Paragraph createTreatmentIcons(@NotNull String allDrugs) {
        String[] drugs = allDrugs.split(Pattern.quote(TREATMENT_DELIMITER));
        Paragraph p = new Paragraph();
        for (String drug : drugs) {
            p.add(Icon.createTreatmentIcon(drug.trim()));
        }
        return p;
    }

    @NotNull
    static Paragraph createChapterFootnoteOffLabel() {
        return new Paragraph().setKeepTogether(true)
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING)
                .add("The Clinical Knowledgebase (CKB) is used to "
                        + "annotate variants of all types with clinical evidence. Only treatment associated evidence with "
                        + "evidence levels \n( ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_A))
                .add(" FDA approved therapy and/or guidelines; ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_B))
                .add(" late clinical trials; ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_C))
                .add(" early clinical trials) can be reported.")
                .add(" Potential evidence items with evidence level \n( ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_D))
                .add(" case reports and preclinical evidence) are not reported.")
                .addStyle(ReportResources.subTextStyle());
    }

    @NotNull
    static Paragraph createChapterFootnoteOnLabel() {
        return new Paragraph().setKeepTogether(true)
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING)
                .add("The Clinical Knowledgebase (CKB) is used to "
                        + "annotate variants of all types with clinical evidence. Only treatment associated evidence with "
                        + "evidence levels \n( ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_A))
                .add(" FDA approved therapy and/or guidelines; ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_B))
                .add(" late clinical trials; ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_C))
                .add(" early clinical trials) can be reported.")
                .add(" Potential evidence items with evidence level \n( ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_D))
                .add(" case reports and preclinical evidence) are not reported.")
                .add("\n")
                .add("\n")
                .add("The iClusion knowledgebase is used to annotate DNA aberrations for potential clinical study eligibility. Of note, "
                        + "clinical study eligibility depends on multiple patient and tumor characteristics of which only the DNA "
                        + "aberrations are considered in this report.")
                .addStyle(ReportResources.subTextStyle());
    }
}
