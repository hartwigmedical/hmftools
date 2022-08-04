package com.hartwig.hmftools.orange.report.chapters;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.KnowledgebaseSource;
import com.hartwig.hmftools.common.protect.EvidenceComparator;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportConfig;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.kernel.pdf.action.PdfAction;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;

import org.apache.commons.compress.utils.Lists;
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

        if (reportConfig.maxEvidenceLevel() != null) {
            String maxEvidenceLevelString = reportConfig.maxEvidenceLevel().toString();
            document.add(note(" * Treatments are reported up to a maximum evidence level of " + maxEvidenceLevelString + "."));
        }

        addTreatmentSection(document, "Applicable", report.protect().reportableEvidences(), report.protect().reportableTrials());

        addTreatmentSection(document,
                "Other potentially interesting",
                report.protect().unreportedEvidences(),
                report.protect().unreportedTrials());
    }

    @NotNull
    private static Paragraph note(@NotNull String message) {
        return new Paragraph(message).addStyle(ReportResources.subTextStyle());
    }

    private void addTreatmentSection(@NotNull Document document, @NotNull String header, @NotNull List<ProtectEvidence> evidences,
            @NotNull List<ProtectEvidence> trials) {
        Map<String, List<ProtectEvidence>> onLabelTreatments = buildEvidenceByEventMap(evidences, true);
        Map<String, List<ProtectEvidence>> offLabelTreatments = buildEvidenceByEventMap(evidences, false);
        document.add(createTreatmentTable(header + " on-label evidence", onLabelTreatments));
        document.add(createTreatmentTable(header + " off-label evidence", offLabelTreatments));

        Map<String, List<ProtectEvidence>> onLabelTrials = buildEvidenceByEventMap(trials, true);
        document.add(createTreatmentTable(header + " trials", onLabelTrials));
    }

    @NotNull
    private static Map<String, List<ProtectEvidence>> buildEvidenceByEventMap(@NotNull List<ProtectEvidence> evidences,
            boolean requireOnLabel) {
        Map<String, List<ProtectEvidence>> evidencePerEventMap = Maps.newHashMap();

        for (ProtectEvidence evidence : evidences) {
            if (evidence.onLabel() == requireOnLabel) {
                String event = evidence.gene() != null ? evidence.gene() + " " + evidence.event() : evidence.event();
                if (event.contains("p.")) {
                    event = AminoAcids.forceSingleLetterProteinAnnotation(event);
                }

                List<ProtectEvidence> eventEvidences = evidencePerEventMap.get(event);
                if (eventEvidences == null) {
                    eventEvidences = Lists.newArrayList();
                }

                if (!containsEvidenceWithDisplay(eventEvidences, evidence)) {
                    eventEvidences.add(evidence);
                }

                evidencePerEventMap.put(event, eventEvidences);
            }
        }
        return evidencePerEventMap;
    }

    private static boolean containsEvidenceWithDisplay(@NotNull List<ProtectEvidence> evidences, @NotNull ProtectEvidence evidenceToMatch) {
        for (ProtectEvidence evidence : evidences) {
            if (display(evidence).equals(display(evidenceToMatch))) {
                return true;
            }
        }
        return false;
    }

    @NotNull
    private Table createTreatmentTable(@NotNull String title, @NotNull Map<String, List<ProtectEvidence>> treatmentByEventMap) {
        Table treatmentTable = Tables.createContent(contentWidth(),
                new float[] { 1, 1, 1 },
                new Cell[] { Cells.createHeader("Event"), Cells.createHeader("Responsive Evidence"),
                        Cells.createHeader("Resistance Evidence") });

        EvidenceLevel maxReportingLevel = reportConfig.maxEvidenceLevel();

        boolean hasEvidence = false;
        for (EvidenceLevel level : EvidenceLevel.values()) {
            if (maxReportingLevel == null || !maxReportingLevel.isHigher(level)) {
                if (addEvidenceWithMaxLevel(treatmentTable, treatmentByEventMap, level)) {
                    hasEvidence = true;
                }
            }
        }

        if (hasEvidence) {
            return Tables.createWrapping(treatmentTable, title);
        } else {
            return Tables.createEmpty(title, contentWidth());
        }
    }

    private boolean addEvidenceWithMaxLevel(@NotNull Table table, @NotNull Map<String, List<ProtectEvidence>> treatmentByEventMap,
            @NotNull EvidenceLevel allowedHighestLevel) {
        Set<String> sortedEvents = Sets.newTreeSet(treatmentByEventMap.keySet());
        boolean hasEvidence = false;
        for (String event : sortedEvents) {
            List<ProtectEvidence> evidences = treatmentByEventMap.get(event);
            if (allowedHighestLevel == highestEvidence(evidences) && containsEvidenceForDisplay(evidences)) {
                table.addCell(Cells.createContent(event));

                Table responsiveTable = Tables.createContent(contentWidth() / 3, new float[] { 1 }, new Cell[] {});
                for (ProtectEvidence responsive : filterOnDirections(evidences, RESPONSIVE_DIRECTIONS)) {
                    Cell cell = Cells.createTransparent(display(responsive));
                    String url = url(responsive);
                    if (!url.isEmpty()) {
                        cell.addStyle(ReportResources.urlStyle()).setAction(PdfAction.createURI(url));
                    }
                    responsiveTable.addCell(cell);
                }
                table.addCell(Cells.createContent(responsiveTable));

                Table resistantTable = Tables.createContent(contentWidth() / 3, new float[] { 1 }, new Cell[] {});
                for (ProtectEvidence resistant : filterOnDirections(evidences, RESISTANT_DIRECTIONS)) {
                    Cell cell = Cells.createTransparent(display(resistant));
                    String url = url(resistant);
                    if (!url.isEmpty()) {
                        cell.addStyle(ReportResources.urlStyle()).setAction(PdfAction.createURI(url));
                    }
                    resistantTable.addCell(cell);
                }
                table.addCell(Cells.createContent(resistantTable));
                hasEvidence = true;
            }
        }

        return hasEvidence;
    }

    private static boolean containsEvidenceForDisplay(@NotNull List<ProtectEvidence> evidences) {
        for (ProtectEvidence evidence : evidences) {
            if (evidence.direction().isResponsive() || evidence.direction().isResistant()) {
                return true;
            }
        }
        return false;
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
        Set<String> sourceNames = Sets.newHashSet();
        for (KnowledgebaseSource source : evidence.sources()) {
            sourceNames.add(source.name().display());
        }

        StringJoiner sources = new StringJoiner(", ");
        for (String sourceName : sourceNames) {
            sources.add(sourceName);
        }
        return evidence.treatment() + " (" + evidence.level() + " - " + sources + ")";
    }

    @NotNull
    private static Set<ProtectEvidence> filterOnDirections(@NotNull List<ProtectEvidence> evidences,
            @NotNull Set<EvidenceDirection> allowedDirections) {
        Set<ProtectEvidence> filtered = Sets.newTreeSet(new EvidenceComparator());
        for (ProtectEvidence evidence : evidences) {
            if (allowedDirections.contains(evidence.direction())) {
                filtered.add(evidence);
            }
        }
        return filtered;
    }

    @NotNull
    private static String url(@NotNull ProtectEvidence evidence) {
        String urlString = Strings.EMPTY;
        for (KnowledgebaseSource source : evidence.sources()) {

            if (source.evidenceUrls().isEmpty()) {
                urlString = Strings.EMPTY;
            }

            for (String url : source.evidenceUrls()) {
                if (url.contains("pubmed")) {
                    // We prefer pubmed URLs over all other URLs so if there is one pubmed then we use that.
                    urlString = url;
                } else if (url.contains("google")) {
                    // If there are no pubmeds, and the first url refers to google we remove it.
                    urlString = Strings.EMPTY;
                } else {
                    urlString = url;
                }
            }

        }
        return urlString;
    }
}