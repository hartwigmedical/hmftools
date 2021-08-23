package com.hartwig.hmftools.orange.report.chapters;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.algo.selection.EvidenceSelector;
import com.hartwig.hmftools.orange.report.ReportConfig;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.TableUtil;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.kernel.pdf.action.PdfAction;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

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

    @NotNull
    @Override
    public PageSize pageSize() {
        return PageSize.A4;
    }

    @Override
    public void render(@NotNull final Document document) {
        document.add(new Paragraph(name()).addStyle(ReportResources.chapterTitleStyle()));

        if (!reportConfig.reportGermline()) {
            document.add(note(" * Evidence from germline events is filtered"));
        }

        if (reportConfig.maxEvidenceLevel() != null) {
            String maxEvidenceLevelString = reportConfig.maxEvidenceLevel().toString();
            document.add(note(" * Treatments are reported up to a maximum evidence level of " + maxEvidenceLevelString + "."));
        }

        List<ProtectEvidence> reported = EvidenceSelector.reported(report.protect());
        addTreatmentSection(document, "Reported", reported);

        List<ProtectEvidence> unreported = EvidenceSelector.unreported(report.protect());
        addTreatmentSection(document, "Unreported", unreported);
    }

    private void addTreatmentSection(@NotNull Document document, @NotNull String header, @NotNull List<ProtectEvidence> evidences) {
        List<ProtectEvidence> noIclusion = EvidenceSelector.noIclusion(evidences);
        Map<String, List<ProtectEvidence>> onLabelTreatments =
                EvidenceSelector.buildTreatmentMap(noIclusion, reportConfig.reportGermline(), true);
        Map<String, List<ProtectEvidence>> offLabelTreatments =
                EvidenceSelector.buildTreatmentMap(noIclusion, reportConfig.reportGermline(), false);
        document.add(createTreatmentTable(header + " on-label evidence", onLabelTreatments));
        document.add(createTreatmentTable(header + " off-label evidence", offLabelTreatments));

        List<ProtectEvidence> iclusion = EvidenceSelector.iclusionOnly(evidences);
        Map<String, List<ProtectEvidence>> trials = EvidenceSelector.buildTreatmentMap(iclusion, reportConfig.reportGermline(), true);
        document.add(createTreatmentTable(header + " trials", trials));
    }

    @NotNull
    private static Paragraph note(@NotNull String message) {
        return new Paragraph(message).addStyle(ReportResources.subTextStyle());
    }

    @NotNull
    private Table createTreatmentTable(@NotNull String title, @NotNull Map<String, List<ProtectEvidence>> treatmentMap) {
        Table treatmentTable = TableUtil.createReportContentTable(contentWidth(),
                new float[] { 1, 1, 1 },
                new Cell[] { TableUtil.createHeaderCell("Treatment"), TableUtil.createHeaderCell("Responsive Evidence"),
                        TableUtil.createHeaderCell("Resistance Evidence") });

        EvidenceLevel maxReportingLevel = reportConfig.maxEvidenceLevel();

        boolean hasEvidence = false;
        for (EvidenceLevel level : EvidenceLevel.values()) {
            if (maxReportingLevel == null || !maxReportingLevel.isHigher(level)) {
                if (addEvidenceWithMaxLevel(treatmentTable, treatmentMap, level)) {
                    hasEvidence = true;
                }
            }
        }

        if (hasEvidence) {
            return TableUtil.createWrappingReportTable(treatmentTable, title);
        } else {
            return TableUtil.createEmptyTable(title, contentWidth());
        }
    }

    private static boolean addEvidenceWithMaxLevel(@NotNull Table table, @NotNull Map<String, List<ProtectEvidence>> treatmentMap,
            @NotNull EvidenceLevel allowedHighestLevel) {
        Set<String> sortedTreatments = Sets.newTreeSet(treatmentMap.keySet());
        boolean hasEvidence = false;
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
