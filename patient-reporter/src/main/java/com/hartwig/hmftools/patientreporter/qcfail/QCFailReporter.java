package com.hartwig.hmftools.patientreporter.qcfail;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Optional;

import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocationFunctions;
import com.hartwig.hmftools.common.lims.Lims;
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
        QCFailStudy study = QCFailStudy.fromSampleId(sampleMetadata.tumorSampleId());

        if (study == null) {
            throw new IllegalStateException("Could not derive study for QC fail report for " + sampleMetadata.tumorSampleId());
        }

        PatientTumorLocation patientTumorLocation =
                PatientTumorLocationFunctions.findPatientTumorLocationForSample(reportData.patientTumorLocations(),
                        sampleMetadata.tumorSampleId());

        String wgsPurityString = null;
        if (reason == QCFailReason.BELOW_DETECTION_THRESHOLD || reason == QCFailReason.POST_ANALYSIS_FAIL) {
            // In these cases we have done full WGS.
            PurityContext purityContext = FittedPurityFile.read(purplePurityTsv);
            boolean hasReliablePurity = CheckPurpleQuality.checkHasReliablePurity(purityContext);

            wgsPurityString =
                    hasReliablePurity ? new DecimalFormat("#'%'").format(purityContext.bestFit().purity() * 100) : Lims.PURITY_NOT_RELIABLE_STRING;

        }

        SampleReport sampleReport = SampleReportFactory.fromLimsModel(sampleMetadata, reportData.limsModel(), patientTumorLocation);

        return ImmutableQCFailReport.builder()
                .sampleReport(sampleReport)
                .reason(reason)
                .study(study)
                .wgsPurityString(wgsPurityString)
                .comments(Optional.ofNullable(comments))
                .isCorrectedReport(correctedReport)
                .signaturePath(reportData.signaturePath())
                .logoRVAPath(reportData.logoRVAPath())
                .logoCompanyPath(reportData.logoCompanyPath())
                .build();
    }
}
