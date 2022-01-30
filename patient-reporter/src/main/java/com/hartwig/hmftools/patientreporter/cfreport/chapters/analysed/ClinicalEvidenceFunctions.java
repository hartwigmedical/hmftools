package com.hartwig.hmftools.patientreporter.cfreport.chapters.analysed;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceType;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.Icon;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.cfreport.data.ClinicalTrials;
import com.hartwig.hmftools.patientreporter.cfreport.data.EvidenceItems;
import com.itextpdf.kernel.pdf.action.PdfAction;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.element.Text;
import com.itextpdf.layout.property.VerticalAlignment;

import org.apache.commons.compress.utils.Lists;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class ClinicalEvidenceFunctions {

    private ClinicalEvidenceFunctions() {
    }

    private static final String TREATMENT_DELIMITER = " + ";

    private static final String RESPONSE_SYMBOL = "\u25B2";
    private static final String RESISTENT_SYMBOL = "\u25BC";
    private static final String PREDICTED_SYMBOL = "P";

    private static final Set<EvidenceDirection> RESISTANT_DIRECTIONS =
            Sets.newHashSet(EvidenceDirection.RESISTANT, EvidenceDirection.PREDICTED_RESISTANT);
    private static final Set<EvidenceDirection> RESPONSE_DIRECTIONS =
            Sets.newHashSet(EvidenceDirection.RESPONSIVE, EvidenceDirection.PREDICTED_RESPONSIVE);
    private static final Set<EvidenceDirection> PREDICTED =
            Sets.newHashSet(EvidenceDirection.PREDICTED_RESISTANT, EvidenceDirection.PREDICTED_RESPONSIVE);

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
            if (evidence.treatment().equals(evidenceToCheck.treatment()) && StringUtils.equals(evidence.gene(), evidenceToCheck.gene())
                    && evidence.event().equals(evidenceToCheck.event())) {
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
                new float[] { 25, 140, 45, 25, 40, 140, 60 },
                new Cell[] { TableUtil.createHeaderCell("Treatment", 2), TableUtil.createHeaderCell("Match", 1),
                        TableUtil.createHeaderCell("Level", 1), TableUtil.createHeaderCell("Response", 1),
                        TableUtil.createHeaderCell("Genomic event", 1), TableUtil.createHeaderCell("Evidence links", 1) });

        treatmentTable = addingDataIntoTable(treatmentTable, treatmentMap, title, contentWidth, "evidence");
        return treatmentTable;
    }

    @NotNull
    public static Table createTrialTable(@NotNull String title, @NotNull Map<String, List<ProtectEvidence>> treatmentMap,
            float contentWidth) {
        Table treatmentTable = TableUtil.createReportContentTable(contentWidth,
                new float[] { 20, 180, 45, 180 },
                new Cell[] { TableUtil.createHeaderCell("Trial", 2), TableUtil.createHeaderCell("Match", 1),
                        TableUtil.createHeaderCell("Genomic event", 1) });

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
            return TableUtil.createWrappingReportTable(title, null, treatmentTable);
        } else {
            return TableUtil.createNoneReportTable(title, null);
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
            @NotNull EvidenceLevel allowedHighestLevel, @NotNull String evidenceType) {
        Set<String> sortedTreatments = Sets.newTreeSet(treatmentMap.keySet().stream().collect(Collectors.toSet()));
        boolean hasEvidence = false;
        for (String treatment : sortedTreatments) {
            List<ProtectEvidence> evidences = treatmentMap.get(treatment);
            if (allowedHighestLevel == highestEvidence(treatmentMap.get(treatment))) {
                table.addCell(TableUtil.createContentCell(createTreatmentIcons(treatment)).setVerticalAlignment(VerticalAlignment.TOP));
                table.addCell(TableUtil.createContentCell(treatment));

                Table typeTable = new Table(new float[] { 1 });
                Table levelTable = new Table(new float[] { 1 });
                Table responseTable = new Table(new float[] { 1, 1 });

                Table responsiveTable = new Table(new float[] { 1 });
                Table linksTable = new Table(new float[] { 1 });

                for (ProtectEvidence responsive : sort(evidences)) {
                    Cell cellGenomic = TableUtil.createTransparentCell(display(responsive));

                    if (!evidenceType.equals("trial")) {
                        cellGenomic.addStyle(ReportResources.urlStyle())
                                .setAction(PdfAction.createURI("https://ckbhome.jax.org/gene/grid"));
                    } else {
                        cellGenomic.addStyle(ReportResources.urlStyle())
                                .setAction(PdfAction.createURI(ClinicalTrials.createLinkiClusion(responsive)));
                    }

                    Cell cellType;
                    cellType = TableUtil.createTransparentCell(new Paragraph(determineEvidenceType(responsive)));
                    typeTable.addCell(cellType);

                    Cell cellLevel;
                    Cell cellPredicted = TableUtil.createTransparentCell(Strings.EMPTY);
                    Cell cellResistent = TableUtil.createTransparentCell(Strings.EMPTY);
                    if (!evidenceType.equals("trial")) {

                        if (PREDICTED.contains(responsive.direction())) {
                            cellPredicted = TableUtil.createTransparentCell(PREDICTED_SYMBOL).addStyle(ReportResources.predictedStyle());
                        }

                        if (RESISTANT_DIRECTIONS.contains(responsive.direction())) {
                            cellResistent = TableUtil.createTransparentCell(RESISTENT_SYMBOL).addStyle(ReportResources.resistentStyle());
                        }

                        if (RESPONSE_DIRECTIONS.contains(responsive.direction())) {
                            cellResistent = TableUtil.createTransparentCell(RESPONSE_SYMBOL).addStyle(ReportResources.responseStyle());
                        }

                        cellLevel = TableUtil.createTransparentCell(new Paragraph(Icon.createLevelIcon(responsive.level().name())));

                        levelTable.addCell(cellLevel);
                        responseTable.addCell(cellResistent);
                        responseTable.addCell(cellPredicted);

                    }

                    responsiveTable.addCell(cellGenomic);

                    Cell publications = TableUtil.createTransparentCell(Strings.EMPTY);
                    if (evidenceType.equals("evidence")) {
                        publications = TableUtil.createTransparentCell(EvidenceItems.createLinksPublications(responsive));
                        linksTable.addCell(publications);
                    } else {
                        linksTable.addCell(publications);
                    }
                }

                if (evidenceType.equals("evidence")) {
                    table.addCell(TableUtil.createContentCell(typeTable));
                    table.addCell(TableUtil.createContentCell(levelTable));
                    table.addCell(TableUtil.createContentCell(responseTable));
                    table.addCell(TableUtil.createContentCell(responsiveTable));
                    table.addCell(TableUtil.createContentCell(linksTable));
                } else {
                    table.addCell(TableUtil.createContentCell(typeTable));
                    table.addCell(TableUtil.createContentCell(responsiveTable));
                }

                hasEvidence = true;
            }
        }

        return hasEvidence;
    }

    @NotNull
    private static String determineEvidenceType(@NotNull ProtectEvidence evidence) {

        String evidenceRank = Strings.EMPTY;
        String evidenceSource = evidence.evidenceType().display();
        if (evidence.evidenceType().equals(ProtectEvidenceType.CODON_MUTATION) || evidence.evidenceType()
                .equals(ProtectEvidenceType.EXON_MUTATION)) {
            evidenceRank = String.valueOf(evidence.rangeRank());
        }

        String evidenceMerged;
        if (!evidenceRank.isEmpty()) {
            evidenceMerged = evidenceSource + " " + evidenceRank;
        } else {
            evidenceMerged = evidenceSource;
        }
        return evidenceMerged;
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
        String event = evidence.gene() != null ? evidence.gene() + " " + evidence.event() : evidence.event();
        if (event.contains("p.")) {
            event = AminoAcids.forceSingleLetterProteinAnnotation(event);
        }

        return event;
    }

    @NotNull
    public static List<ProtectEvidence> sort(@NotNull List<ProtectEvidence> evidenceItems) {
        return evidenceItems.stream().sorted((item1, item2) -> {
            if (item1.treatment().equals(item2.treatment())) {
                if (item1.level().equals(item2.level())) {
                    if (item1.direction().equals(item2.direction())) {
                        return item1.direction().compareTo(item2.direction());
                    } else {
                        return item1.direction().compareTo(item2.direction());
                    }
                } else {
                    return item1.level().compareTo(item2.level());
                }
            } else {
                return item1.treatment().compareTo(item2.treatment());
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    public static Paragraph noteGlossaryTerms() {
        return new Paragraph("The symbol ( ").add(new Text(RESPONSE_SYMBOL).addStyle(ReportResources.responseStyle()))
                .add(" ) means that the evidence is responsive. The symbol ( ")
                .add(new Text(RESISTENT_SYMBOL).addStyle(ReportResources.resistentStyle()))
                .add(" ) means that the evidence is resistant. The abbreviation ( ")
                .add(new Text(PREDICTED_SYMBOL).addStyle(ReportResources.predictedStyle()))
                .add(" mentioned after the level of evidence) indicates the evidence is predicted "
                        + "responsive/resistent. More details about CKB can be found in their")
                .addStyle(ReportResources.subTextStyle())
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING)
                .add(new Text(" Glossary Of Terms").addStyle(ReportResources.urlStyle())
                        .setAction(PdfAction.createURI("https://ckbhome.jax.org/about/glossaryOfTerms")))
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING);
    }

    @NotNull
    public static Paragraph noteEvidence() {
        return new Paragraph().setFixedLeading(ReportResources.BODY_TEXT_LEADING)
                .add("The Clinical Knowledgebase (CKB) is used to annotate variants of all types with clinical evidence. "
                        + "Only treatment associated evidence with evidence levels ( \n( ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_A))
                .add(" FDA approved therapy and/or guidelines; ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_B))
                .add(" late clinical trials; ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_C))
                .add(" early clinical trials) can be reported. Potential evidence items with evidence level  \n( ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_D))
                .add(" case reports and preclinical evidence) are not reported.")
                .addStyle(ReportResources.subTextStyle());
    }

    @NotNull
    public static Paragraph note(@NotNull String message) {
        return new Paragraph(message).addStyle(ReportResources.subTextStyle());
    }
}
