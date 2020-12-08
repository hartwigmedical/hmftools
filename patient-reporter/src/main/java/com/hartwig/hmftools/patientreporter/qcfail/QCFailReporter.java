package com.hartwig.hmftools.patientreporter.qcfail;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Optional;

import com.hartwig.hmftools.common.clinical.PatientPrimaryTumor;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumorFunctions;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsCohort;
import com.hartwig.hmftools.common.lims.LimsStudy;
import com.hartwig.hmftools.common.purple.CheckPurpleQuality;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContextFile;
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
            @NotNull String purpleQCFile, @Nullable String comments, boolean correctedReport) throws IOException {

        PatientPrimaryTumor patientPrimaryTumor =
                PatientPrimaryTumorFunctions.findPrimaryTumorForSample(reportData.patientPrimaryTumors(),
                        sampleMetadata.tumorSampleId());
        SampleReport sampleReport = SampleReportFactory.fromLimsModel(sampleMetadata, reportData.limsModel(), patientPrimaryTumor);

        LimsCohort cohort = sampleReport.cohort();

        if (cohort == cohort.NON_CANCER) {
            throw new IllegalStateException("QC fail report not supported for non-cancer study samples: " + sampleMetadata.tumorSampleId());
        }

        String wgsPurityString = null;
        if (reason.isDeepWGSDataAvailable()) {
            PurityContext purityContext = PurityContextFile.readWithQC(purpleQCFile, purplePurityTsv);
            String formattedPurity = new DecimalFormat("#'%'").format(purityContext.bestFit().purity() * 100);
            boolean hasReliablePurity = CheckPurpleQuality.checkHasReliablePurity(purityContext);

            wgsPurityString = hasReliablePurity ? formattedPurity : Lims.PURITY_NOT_RELIABLE_STRING;
        }


        return ImmutableQCFailReport.builder()
                .sampleReport(sampleReport)
                .qsFormNumber(reason.qcFormNumber())
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
