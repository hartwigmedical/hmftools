package com.hartwig.hmftools.orange.report.chapters;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.AminoAcidFunctions;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportConfig;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.DocumentUtil;
import com.hartwig.hmftools.orange.report.util.TableUtil;
import com.itextpdf.kernel.pdf.action.PdfAction;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ClinicalEvidenceChapter implements ReportChapter {

    private static final Set<EvidenceDirection> RESPONSIVE_DIRECTIONS =
            Sets.newHashSet(EvidenceDirection.RESPONSIVE, EvidenceDirection.PREDICTED_RESPONSIVE);
    private static final Set<EvidenceDirection> RESISTANT_DIRECTIONS =
            Sets.newHashSet(EvidenceDirection.RESISTANT, EvidenceDirection.PREDICTED_RESISTANT);

    @NotNull
    private final OrangeReport report;
    @NotNull
    private final ReportConfig reportConfig;

    public ClinicalEvidenceChapter(@NotNull final OrangeReport report, @NotNull final ReportConfig reportConfig) {
        this.report = report;
        this.reportConfig = reportConfig;
    }

    @NotNull
    @Override
    public String name() {
        return "Clinical Evidence";
    }

    @Override
    public void render(@NotNull final Document document) {
        document.add(new Paragraph("Clinical Evidence").addStyle(ReportResources.chapterTitleStyle()));

        if (!reportConfig.reportGermline()) {
            document.add(note(" * Evidence from germline events is filtered"));
        }

        if (reportConfig.maxReportingLevel() != null) {
            document.add(note(" * Treatments are reported up to a maximum evidence level of '" + reportConfig.maxReportingLevel().toString()
                    + "'.'"));
        }

        List<ProtectEvidence> noIclusion = noIclusion(report.protect());

        String onLabelTitle = "On-Label Evidence";
        Table onLabelTable = createTreatmentTable(onLabelTitle, toTreatmentMap(noIclusion, reportConfig.reportGermline(), true));

        String offLabelTitle = "Off-Label Evidence";
        Table offLabelTable = createTreatmentTable(offLabelTitle, toTreatmentMap(noIclusion, reportConfig.reportGermline(), false));

        List<ProtectEvidence> iclusion = iclusionOnly(report.protect());
        String trialTitle = "Trials";
        Table trialTable = createTreatmentTable(trialTitle, toTreatmentMap(iclusion, reportConfig.reportGermline(), true));

        DocumentUtil.addCheckedTable(document, onLabelTitle, onLabelTable);
        DocumentUtil.addCheckedTable(document, offLabelTitle, offLabelTable);
        DocumentUtil.addCheckedTable(document, trialTitle, trialTable);
    }

    @NotNull
    private static List<ProtectEvidence> iclusionOnly(@NotNull List<ProtectEvidence> evidences) {
        List<ProtectEvidence> filtered = Lists.newArrayList();
        for (ProtectEvidence evidence : evidences) {
            if (evidence.sources().contains(Knowledgebase.ICLUSION)) {
                filtered.add(evidence);
            }
        }
        return filtered;
    }

    @NotNull
    private static List<ProtectEvidence> noIclusion(@NotNull List<ProtectEvidence> evidences) {
        List<ProtectEvidence> filtered = Lists.newArrayList();
        for (ProtectEvidence evidence : evidences) {
            if (!evidence.sources().contains(Knowledgebase.ICLUSION)) {
                filtered.add(evidence);
            }
        }
        return filtered;
    }

    @NotNull
    private static Paragraph note(@NotNull String message) {
        return new Paragraph(message).addStyle(ReportResources.subTextStyle());
    }

    @NotNull
    private static Map<String, List<ProtectEvidence>> toTreatmentMap(@NotNull List<ProtectEvidence> evidences, boolean reportGermline,
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

    @Nullable
    private Table createTreatmentTable(@NotNull String title, @NotNull Map<String, List<ProtectEvidence>> treatmentMap) {
        if (treatmentMap.isEmpty()) {
            return null;
        }

        Table treatmentTable = TableUtil.createReportContentTable(new float[] { 1, 1, 1 },
                new Cell[] { TableUtil.createHeaderCell("Treatment"), TableUtil.createHeaderCell("Responsive Evidence"),
                        TableUtil.createHeaderCell("Resistance Evidence") });

        EvidenceLevel maxReportingLevel = reportConfig.maxReportingLevel();

        for (EvidenceLevel level : EvidenceLevel.values()) {
            if (maxReportingLevel == null || !maxReportingLevel.isHigher(level)) {
                addEvidenceWithMaxLevel(treatmentTable, treatmentMap, level);
            }
        }

        return TableUtil.createWrappingReportTable(treatmentTable, title);
    }

    private static void addEvidenceWithMaxLevel(@NotNull Table table, @NotNull Map<String, List<ProtectEvidence>> treatmentMap,
            @NotNull EvidenceLevel allowedHighestLevel) {
        Set<String> sortedTreatments = Sets.newTreeSet(treatmentMap.keySet());
        for (String treatment : sortedTreatments) {
            List<ProtectEvidence> evidences = treatmentMap.get(treatment);
            if (allowedHighestLevel == highestEvidence(treatmentMap.get(treatment))) {
                table.addCell(TableUtil.createContentCell(treatment));

                Table responsiveTable = new Table(1);
                for (ProtectEvidence responsive : filterOnDirections(evidences, RESPONSIVE_DIRECTIONS)) {
                    Cell cell = TableUtil.createTransparentCell(display(responsive));
                    String url = url(responsive);
                    if (!url.isEmpty()) {
                        cell.addStyle(ReportResources.urlStyle()).setAction(PdfAction.createURI(url));
                    }
                    responsiveTable.addCell(cell);
                }
                table.addCell(TableUtil.createContentCell(responsiveTable));

                Table resistantTable = new Table(1);
                for (ProtectEvidence resistant : filterOnDirections(evidences, RESISTANT_DIRECTIONS)) {
                    Cell cell = TableUtil.createTransparentCell(display(resistant));
                    String url = url(resistant);
                    if (!url.isEmpty()) {
                        cell.addStyle(ReportResources.urlStyle()).setAction(PdfAction.createURI(url));
                    }
                    resistantTable.addCell(cell);
                }
                table.addCell(TableUtil.createContentCell(resistantTable));
            }
        }
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
            event = AminoAcidFunctions.forceSingleLetterProteinAnnotation(event);
        }
        StringJoiner sources = new StringJoiner(", ");
        for (Knowledgebase source : evidence.sources()) {
            sources.add(source.reportDisplay());
        }
        return event + " (" + evidence.level().toString() + " - " + sources.toString() + ")";
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
}
