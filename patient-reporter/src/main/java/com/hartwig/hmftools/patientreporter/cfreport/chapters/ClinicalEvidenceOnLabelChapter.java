package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.patientreporter.algo.GenomicAnalysis;
import com.itextpdf.layout.Document;

import org.jetbrains.annotations.NotNull;

public class ClinicalEvidenceOnLabelChapter implements ReportChapter {

    @Override
    @NotNull
    public String name() {
        return "Therapy details (Tumor type specific)";
    }

    @NotNull
    private final GenomicAnalysis analysis;

    public ClinicalEvidenceOnLabelChapter(@NotNull final GenomicAnalysis analysis) {
        this.analysis = analysis;
    }


    @Override
    public void render(@NotNull final Document document) {

        List<ProtectEvidence> reportedOnLabel = analysis.tumorSpecificEvidence();
        addTreatmentSection(document, "Tumor type specific evidence", reportedOnLabel);

        List<ProtectEvidence> reportedStudies = analysis.clinicalTrials();
        addTrialSection(document, "Tumor type specific clinical trials (NL)", reportedStudies);
        document.add(ClinicalEvidenceFunctions.note("Potential eligibility for DRUP is dependent on tumor type details "
                + "therefore certain tumor types may not be eligible for the DRUP. "));
        document.add(ClinicalEvidenceFunctions.note(""));
        document.add(ClinicalEvidenceFunctions.note("The iClusion knowledgebase is used to annotate DNA aberrations for potential "
                + "clinical study eligibility. Of note, clinical study eligibility depends on multiple patient and tumor "
                + "characteristics of which only the DNA aberrations are considered in this report."));
        document.add(ClinicalEvidenceFunctions.note(""));
        document.add(ClinicalEvidenceFunctions.noteEvidence());
    }

    private void addTreatmentSection(@NotNull Document document, @NotNull String header, @NotNull List<ProtectEvidence> evidences) {
        Map<String, List<ProtectEvidence>> onLabelTreatments =
                ClinicalEvidenceFunctions.buildTreatmentMap(evidences, true, true);
        document.add(ClinicalEvidenceFunctions.createTreatmentTable(header, onLabelTreatments, contentWidth()));
    }

    private void addTrialSection(@NotNull Document document, @NotNull String header, @NotNull List<ProtectEvidence> evidences) {
        Map<String, List<ProtectEvidence>> onLabelTreatments =
                ClinicalEvidenceFunctions.buildTreatmentMap(evidences, true, true);
        document.add(ClinicalEvidenceFunctions.createTrialTable(header, onLabelTreatments, contentWidth()));
    }
}