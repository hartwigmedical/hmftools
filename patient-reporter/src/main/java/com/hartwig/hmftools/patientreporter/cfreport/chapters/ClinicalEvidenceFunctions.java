package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import java.net.URL;
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

    private static final Set<EvidenceDirection> RESISTANT_DIRECTIONS =
            Sets.newHashSet(EvidenceDirection.RESISTANT, EvidenceDirection.PREDICTED_RESISTANT);
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
                new float[] { 25, 160, 25, 35, 160 },
                new Cell[] { TableUtil.createHeaderCell("Treatment", 2), TableUtil.createHeaderCell("Level", 1),
                        TableUtil.createHeaderCell("Response", 1), TableUtil.createHeaderCell("Genomic event", 1) });

        treatmentTable = addingDataIntoTable(treatmentTable, treatmentMap, title, contentWidth, "evidence");
        return treatmentTable;
    }

    @NotNull
    public static Table createTrialTable(@NotNull String title, @NotNull Map<String, List<ProtectEvidence>> treatmentMap,
            float contentWidth) {
        Table treatmentTable = TableUtil.createReportContentTable(contentWidth,
                new float[] { 20, 180, 180 },
                new Cell[] { TableUtil.createHeaderCell("Trial", 2), TableUtil.createHeaderCell("Genomic event", 1) });

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

                Table levelTable = new Table(new float[] { 1 });
                Table responseTable = new Table(new float[] { 1, 1 });

                Table responsiveTable = new Table(new float[] { 1, 1 });

                for (ProtectEvidence responsive : evidences) {
                    Cell cellGenomic = TableUtil.createTransparentCell(display(responsive));

                    if (!evidenType.equals("trial")) {
                        cellGenomic.addStyle(ReportResources.urlStyle())
                                .setAction(PdfAction.createURI("https://ckbhome.jax.org/gene/grid"));
                    } else {
                        cellGenomic.addStyle(ReportResources.urlStyle())
                                .setAction(PdfAction.createURI(createLinkiClusion(responsive)));
                    }


                    Cell cellLevel;
                    Cell cellPredicted;
                    Cell cellResistent;
                    if (!evidenType.equals("trial")) {
                        String predicted = " ";
                        if (PREDICTED.contains(responsive.direction())) {
                            predicted = "D";
                        }

                        String resistent = " ";
                        if (RESISTANT_DIRECTIONS.contains(responsive.direction())) {
                            resistent = "C";
                        }

                        cellLevel = TableUtil.createTransparentCell(new Paragraph(Icon.createLevelIcon(responsive.level().name())));
                        cellPredicted = TableUtil.createTransparentCell(new Paragraph(Icon.createLevelIcon(predicted)));
                        cellResistent = TableUtil.createTransparentCell(new Paragraph(Icon.createLevelIcon(resistent)));
                        levelTable.addCell(cellLevel);
                        responseTable.addCell(cellResistent);
                        responseTable.addCell(cellPredicted);

                    }

                    responsiveTable.addCell(cellGenomic);

                    if (evidenType.equals("evidence")) {
                        Cell cellnumbers = createLinksPublications(responsive);
                        responsiveTable.addCell(cellnumbers);
                    } else {
                        Cell cellnumbers = TableUtil.createTransparentCell(Strings.EMPTY);
                        responsiveTable.addCell(cellnumbers);
                    }
                }

                if (evidenType.equals("evidence")) {
                    table.addCell(TableUtil.createContentCell(levelTable));
                    table.addCell(TableUtil.createContentCell(responseTable));
                }
                table.addCell(TableUtil.createContentCell(responsiveTable));

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
    private static String createLinkiClusion(@NotNull ProtectEvidence evidence) {
        String link = Strings.EMPTY;
        for (String url : evidence.urls()) {
            if (url.contains("iclusion")) {
                link = url;
            }
        }
        //We assume iClusion has one link
        return link;
    }

    @NotNull
    private static Cell createLinksPublications(@NotNull ProtectEvidence evidence) {
        List<String> urlForSymlinks = Lists.newArrayList();
        for (String url : evidence.urls()) {
            if (!url.contains("google")) {
                urlForSymlinks.add(url);
            }
        }

        Cell cellnumbers = TableUtil.createTransparentCell(Strings.EMPTY);
        String numbers = urlForSymlinks.size() > 0 ? "[" : Strings.EMPTY;
        for (int i = 0; i < urlForSymlinks.size(); i++) {
            numbers = numbers.concat(Integer.toString(i + 1)).concat(", ");
        }
        numbers = numbers.length() > 1 ? numbers.substring(0, numbers.length() - 2).concat("]") : numbers;
        cellnumbers = TableUtil.createTransparentCell(numbers);

        for (int i = 0; i < urlForSymlinks.size(); i++) {
            cellnumbers.addStyle(ReportResources.urlStyle()).setAction(PdfAction.createURI(urlForSymlinks.get(i)));
        }
        return cellnumbers;
    }

    @NotNull
    public static Paragraph noteEvidence() {
        return new Paragraph().setKeepTogether(true)
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING)
                .add("if PRD mentioned after level then the evidence is predicted response/ressitent otherwise it is real "
                        + "responsive/resistent\n\n")
                .add("The Clinical Knowledgebase (CKB) is used to annotate variants of all types with clinical evidence. "
                        + "Only treatment associated evidence with evidence levels ( \n( ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_A))
                .add(" FDA approved therapy and/or guidelines; ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_B))
                .add(" late clinical trials; ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_C))
                .add(" arly clinical trials) can be reported. Potential evidence items with evidence level  \n( ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_D))
                .add(" case reports and preclinical evidence) are not reported.\n\n")
                .addStyle(ReportResources.subTextStyle());
    }

    @NotNull
    public static Paragraph note(@NotNull String message) {
        return new Paragraph(message).addStyle(ReportResources.subTextStyle());
    }
}
