package com.hartwig.hmftools.patientreporter.qcfail;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Optional;

import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocationFunctions;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsStudy;
import com.hartwig.hmftools.common.purple.CheckPurpleQuality;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.patientreporter.SampleMetadata;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.SampleReportFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class QCFailReporter {

    @NotNull
    private final QCFailReportData reportData;

    public QCFailReporter(@NotNull final QCFailReportData reportData) {
        this.reportData = reportData;
    }

    @NotNull
    public QCFailReport run(@NotNull QCFailReason reason, @NotNull SampleMetadata sampleMetadata, @NotNull String purplePurityTsv,
            @Nullable String comments, boolean correctedReport) throws IOException {
        LimsStudy study = LimsStudy.fromSampleId(sampleMetadata.tumorSampleId());

        if (study == LimsStudy.NON_CANCER_STUDY) {
            throw new IllegalStateException("QC fail report not supported for non-cancer study samples: " + sampleMetadata.tumorSampleId());
        }

        PatientTumorLocation patientTumorLocation =
                PatientTumorLocationFunctions.findPatientTumorLocationForSample(reportData.patientTumorLocations(),
                        sampleMetadata.tumorSampleId());

        String wgsPurityString = null;
        if (reason.isDeepWGSDataAvailable()) {
            PurityContext purityContext = FittedPurityFile.read(purplePurityTsv);
            String formattedPurity = new DecimalFormat("#'%'").format(purityContext.bestFit().purity() * 100);
            boolean hasReliablePurity = CheckPurpleQuality.checkHasReliablePurity(purityContext);

            wgsPurityString = hasReliablePurity ? formattedPurity : Lims.PURITY_NOT_RELIABLE_STRING;
        }

        SampleReport sampleReport = SampleReportFactory.fromLimsModel(sampleMetadata, reportData.limsModel(), patientTumorLocation);

        return ImmutableQCFailReport.builder()
                .sampleReport(sampleReport)
                .reason(reason)
                .wgsPurityString(wgsPurityString)
                .comments(Optional.ofNullable(comments))
                .isCorrectedReport(correctedReport)
                .signaturePath(reportData.signaturePath())
                .logoRVAPath(reportData.logoRVAPath())
                .logoCompanyPath(reportData.logoCompanyPath())
                .build();
    }
}
