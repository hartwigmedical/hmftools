package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.patientreporter.algo.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.algo.GenomicAnalysis;
import com.itextpdf.layout.Document;

import org.jetbrains.annotations.NotNull;

public class ClinicalEvidenceOffLabelChapter implements ReportChapter {

    @Override
    @NotNull
    public String name() {
        return "Therapy details (Other tumor types)";
    }

    @NotNull
    private final AnalysedPatientReport report;

    public ClinicalEvidenceOffLabelChapter(@NotNull final AnalysedPatientReport report) {
        this.report = report;
    }

    @Override
    public void render(@NotNull final Document document) {
        GenomicAnalysis analysis = report.genomicAnalysis();
        List<ProtectEvidence> reportedOffLabel = analysis.offLabelEvidence();
        addTreatmentSection(document, "Evidence on other tumor types", reportedOffLabel);
        document.add(ClinicalEvidenceFunctions.noteEvidence());
    }

    private void addTreatmentSection(@NotNull Document document, @NotNull String header, @NotNull List<ProtectEvidence> evidences) {
        boolean reportGermline = report.sampleReport().germlineReportingLevel().equals(LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION);
        boolean requireOnLabel = true;
        Map<String, List<ProtectEvidence>> offLabelTreatments =
                ClinicalEvidenceFunctions.buildTreatmentMap(evidences, reportGermline, requireOnLabel);
        document.add(ClinicalEvidenceFunctions.createTreatmentTable(header, offLabelTreatments, contentWidth()));
    }
}
