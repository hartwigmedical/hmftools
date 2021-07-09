package com.hartwig.hmftools.orange.report.chapters;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.components.TableUtil;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public class ClinicalEvidenceChapter implements ReportChapter {

    @NotNull
    private final OrangeReport report;
    private final boolean reportGermline;

    public ClinicalEvidenceChapter(@NotNull final OrangeReport report, final boolean reportGermline) {
        this.report = report;
        this.reportGermline = reportGermline;
    }

    @NotNull
    @Override
    public String name() {
        return "Clinical Evidence";
    }

    @Override
    public void render(@NotNull final Document document) {
        document.add(new Paragraph("Clinical Evidence").addStyle(ReportResources.chapterTitleStyle()));

        Map<String, List<ProtectEvidence>> onLabelTreatmentMap = toTreatmentMap(report.protect(), reportGermline, true);

        Table treatmentTable = TableUtil.createReportContentTable(new float[] { 1, 1, 1 },
                new Cell[] { TableUtil.createHeaderCell("Treatment"), TableUtil.createHeaderCell("Responsive Evidence"),
                        TableUtil.createHeaderCell("Resistance Evidence") });

        Set<String> sortedTreatments = Sets.newTreeSet(onLabelTreatmentMap.keySet());
        for (String treatment : sortedTreatments) {
            List<ProtectEvidence> evidences = onLabelTreatmentMap.get(treatment);
            treatmentTable.addCell(TableUtil.createContentCell(treatment));
            treatmentTable.addCell(TableUtil.createContentCell(evidences.get(0).genomicEvent()));
            treatmentTable.addCell(TableUtil.createContentCell(evidences.get(0).genomicEvent()));
        }

        document.add(TableUtil.createWrappingReportTable(treatmentTable, "On-Label Evidence"));
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
                if (!hasHigherOrEqualEvidenceForEvent(treatmentEvidences, evidence)) {
                    treatmentEvidences.add(evidence);
                }
                evidencePerTreatmentMap.put(evidence.treatment(), treatmentEvidences);
            }
        }
        return evidencePerTreatmentMap;
    }

    private static boolean hasHigherOrEqualEvidenceForEvent(@NotNull List<ProtectEvidence> evidences,
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
}
