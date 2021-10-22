package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.regex.Pattern;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.Icon;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.itextpdf.kernel.pdf.action.PdfAction;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.VerticalAlignment;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class ClinicalEvidenceFunctions {

    private ClinicalEvidenceFunctions() {
    }

    private static final String TREATMENT_DELIMITER = " + ";

    private static final Set<EvidenceDirection> RESPONSIVE_DIRECTIONS =
            Sets.newHashSet(EvidenceDirection.RESPONSIVE, EvidenceDirection.PREDICTED_RESPONSIVE);
    private static final Set<EvidenceDirection> RESISTANT_DIRECTIONS =
            Sets.newHashSet(EvidenceDirection.RESISTANT, EvidenceDirection.PREDICTED_RESISTANT);

    @NotNull
    public static Map<String, List<ProtectEvidence>> buildTreatmentMap(@NotNull List<ProtectEvidence> evidences, boolean reportGermline,
            boolean requireOnLabel) {
        Map<String, List<ProtectEvidence>> evidencePerTreatmentMap = Maps.newHashMap();

        for (ProtectEvidence evidence : evidences) {
            if ((reportGermline || !evidence.germline()) && evidence.onLabel() == requireOnLabel) {
                List<ProtectEvidence> treatmentEvidences = evidencePerTreatmentMap.get(evidence.treatment());
                if (treatmentEvidences == null) {
                    treatmentEvidences = Lists.newArrayList();
                }
                if (!hasHigherOrEqualEvidenceForEventAndTreatment(treatmentEvidences, evidence)) {
                    treatmentEvidences.add(evidence);
                }
                evidencePerTreatmentMap.put(evidence.treatment(), treatmentEvidences);
            }
        }
        return evidencePerTreatmentMap;
    }

    private static boolean hasHigherOrEqualEvidenceForEventAndTreatment(@NotNull List<ProtectEvidence> evidences,
            @NotNull ProtectEvidence evidenceToCheck) {
        for (ProtectEvidence evidence : evidences) {
            if (evidence.treatment().equals(evidenceToCheck.treatment()) && evidence.genomicEvent()
                    .equals(evidenceToCheck.genomicEvent())) {
                if (!evidenceToCheck.level().isHigher(evidence.level())) {
                    return true;
                }
            }
        }
        return false;
    }

    @NotNull
    public static Table createTreatmentTable(@NotNull String title, @NotNull Map<String, List<ProtectEvidence>> treatmentMap,
            float contentWidth) {
        Table treatmentTable = TableUtil.createReportContentTable(contentWidth,
                new float[] { 25, 150, 150, 150 },
                new Cell[] { TableUtil.createHeaderCell("Treatment", 2), TableUtil.createHeaderCell("Responsive Evidence"),
                        TableUtil.createHeaderCell("Resistance Evidence") });

        treatmentTable = addingDataIntoTable(treatmentTable, treatmentMap, title, contentWidth, "evidence");
        return treatmentTable;
    }

    @NotNull
    public static Table createTrialTable(@NotNull String title, @NotNull Map<String, List<ProtectEvidence>> treatmentMap,
            float contentWidth) {
        Table treatmentTable = TableUtil.createReportContentTable(contentWidth,
                new float[] { 20, 180, 180 },
                new Cell[] { TableUtil.createHeaderCell("Trial", 2), TableUtil.createHeaderCell("Eligibility", 2) });

        treatmentTable = addingDataIntoTable(treatmentTable, treatmentMap, title, contentWidth, "trial");
        return treatmentTable;
    }

    @NotNull
    private static Table addingDataIntoTable(@NotNull Table treatmentTable, @NotNull Map<String, List<ProtectEvidence>> treatmentMap,
            @NotNull String title, float contentWidth, @NotNull String evidenType) {
        boolean hasEvidence = false;
        for (EvidenceLevel level : EvidenceLevel.values()) {
            if (addEvidenceWithMaxLevel(treatmentTable, treatmentMap, level, evidenType)) {
                hasEvidence = true;
            }

        }

        if (hasEvidence) {
            return TableUtil.createWrappingReportTable(title, treatmentTable);
        } else {
            return TableUtil.createEmptyTable(title, contentWidth);
        }
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

    private static boolean addEvidenceWithMaxLevel(@NotNull Table table, @NotNull Map<String, List<ProtectEvidence>> treatmentMap,
            @NotNull EvidenceLevel allowedHighestLevel, @NotNull String evidenType) {
        Set<String> sortedTreatments = Sets.newTreeSet(treatmentMap.keySet());
        boolean hasEvidence = false;
        for (String treatment : sortedTreatments) {
            List<ProtectEvidence> evidences = treatmentMap.get(treatment);
            if (allowedHighestLevel == highestEvidence(treatmentMap.get(treatment))) {
                table.addCell(TableUtil.createContentCell(createTreatmentIcons(treatment)).setVerticalAlignment(VerticalAlignment.TOP));
                table.addCell(TableUtil.createContentCell(treatment));

                Table responsiveTable = new Table(new float[] { 1, 1 });
                for (ProtectEvidence responsive : filterOnDirections(evidences, RESPONSIVE_DIRECTIONS)) {
                    Cell cell = TableUtil.createTransparentCell(display(responsive));
                    String url = url(responsive);
                    if (!url.isEmpty()) {
                        cell.addStyle(ReportResources.urlStyle()).setAction(PdfAction.createURI(url));
                    }
                    responsiveTable.addCell(TableUtil.createContentCell(new Paragraph(Icon.createLevelIcon(responsive.level().name()))));
                    responsiveTable.addCell(cell);
                }
                table.addCell(TableUtil.createContentCell(responsiveTable));

                if (evidenType.equals("evidence")) {
                    Table resistantTable = new Table(new float[] { 1, 1 });
                    for (ProtectEvidence resistant : filterOnDirections(evidences, RESISTANT_DIRECTIONS)) {
                        Cell cell = TableUtil.createTransparentCell(display(resistant));
                        String url = url(resistant);
                        if (!url.isEmpty()) {
                            cell.addStyle(ReportResources.urlStyle()).setAction(PdfAction.createURI(url));
                        }
                        resistantTable.addCell(TableUtil.createContentCell(new Paragraph(Icon.createLevelIcon(resistant.level().name()))));
                        resistantTable.addCell(cell);
                    }
                    table.addCell(TableUtil.createContentCell(resistantTable));
                }

                hasEvidence = true;
            }
        }

        return hasEvidence;
    }

    @NotNull
    private static EvidenceLevel highestEvidence(@NotNull List<ProtectEvidence> evidences) {
        EvidenceLevel highest = null;
        for (ProtectEvidence evidence : evidences) {
            if (highest == null || evidence.level().isHigher(highest)) {
                highest = evidence.level();
            }
        }

        return highest;
    }

    @NotNull
    private static String display(@NotNull ProtectEvidence evidence) {
        String event = evidence.genomicEvent();
        if (event.contains("p.")) {
            event = AminoAcids.forceSingleLetterProteinAnnotation(event);
        }
        StringJoiner sources = new StringJoiner(", ");
        for (Knowledgebase source : evidence.sources()) {
            sources.add(source.reportDisplay());
        }

        return event + " (" + sources + ")";
    }

    @NotNull
    private static Set<ProtectEvidence> filterOnDirections(@NotNull List<ProtectEvidence> evidences,
            @NotNull Set<EvidenceDirection> allowedDirections) {
        Set<ProtectEvidence> filtered = Sets.newTreeSet();
        for (ProtectEvidence evidence : evidences) {
            if (allowedDirections.contains(evidence.direction())) {
                filtered.add(evidence);
            }
        }
        return filtered;
    }

    @NotNull
    private static String url(@NotNull ProtectEvidence evidence) {
        if (evidence.urls().isEmpty()) {
            return Strings.EMPTY;
        }

        // We prefer pubmed URLs over all other URLs so if there is one pubmed then we use that.
        for (String url : evidence.urls()) {
            if (url.contains("pubmed")) {
                return url;
            }
        }

        // If there are no pubmeds, and the first url refers to google we remove it.
        String url = evidence.urls().iterator().next();
        if (url.contains("google")) {
            return Strings.EMPTY;
        } else {
            return url;
        }
    }

    @NotNull
    public static Paragraph noteEvidence() {
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
    public static Paragraph note(@NotNull String message) {
        return new Paragraph(message).addStyle(ReportResources.subTextStyle());
    }
}
