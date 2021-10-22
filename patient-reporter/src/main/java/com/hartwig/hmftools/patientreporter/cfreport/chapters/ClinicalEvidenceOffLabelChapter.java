package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.protect.ProtectEvidence;
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
    private final GenomicAnalysis analysis;

    public ClinicalEvidenceOffLabelChapter(@NotNull final GenomicAnalysis analysis) {
        this.analysis = analysis;
    }

    @Override
    public void render(@NotNull final Document document) {

        List<ProtectEvidence> reportedOffLabel = analysis.offLabelEvidence();
        addTreatmentSection(document, "Evidence on other tumor types", reportedOffLabel);
        document.add(ClinicalEvidenceFunctions.noteEvidence());
    }

    private void addTreatmentSection(@NotNull Document document, @NotNull String header, @NotNull List<ProtectEvidence> evidences) {
        Map<String, List<ProtectEvidence>> offLabelTreatments = ClinicalEvidenceFunctions.buildTreatmentMap(evidences, true, false);
        document.add(ClinicalEvidenceFunctions.createTreatmentTable(header, offLabelTreatments, contentWidth()));
    }
}
