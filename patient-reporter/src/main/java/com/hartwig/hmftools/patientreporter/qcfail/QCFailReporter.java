package com.hartwig.hmftools.patientreporter.qcfail;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Optional;

import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocationFunctions;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.patientreporter.SampleMetadata;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.SampleReportFactory;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class QCFailReporter {

    @NotNull
    private final QCFailReportData reportData;

    public QCFailReporter(@NotNull final QCFailReportData reportData) {
        this.reportData = reportData;
    }

    @NotNull
    public QCFailReport run(@NotNull SampleMetadata sampleMetadata, @NotNull QCFailReason reason, @Nullable String comments,
            boolean correctedReport, @Nullable String purplePurityTsv) throws IOException {
        QCFailStudy study = QCFailStudy.fromSampleId(sampleMetadata.tumorSampleId());

        String purity = Strings.EMPTY;
        if (purplePurityTsv != null) {
            PurityContext purityContext = FittedPurityFile.read(purplePurityTsv);
            purity = new DecimalFormat("#'%'").format(purityContext.bestFit().purity() * 100);
        }

        assert study != null;

        PatientTumorLocation patientTumorLocation =
                PatientTumorLocationFunctions.findPatientTumorLocationForSample(reportData.patientTumorLocations(),
                        sampleMetadata.tumorSampleId());

        String calculatedPurity = Strings.EMPTY;
        if (reason == QCFailReason.BELOW_DETECTION_THRESHOLD || reason == QCFailReason.POST_ANALYSIS_FAIL) {

            //check for sample are present in shallow seq db when executed
            reportData.limsModel().purityShallowSeq(sampleMetadata.tumorSampleBarcode());

            calculatedPurity = purity;
        } else {
            calculatedPurity = reportData.limsModel().purityShallowSeq(sampleMetadata.tumorSampleBarcode());
        }

        SampleReport sampleReport = SampleReportFactory.fromLimsModel(sampleMetadata,
                reportData.limsModel(),
                patientTumorLocation, calculatedPurity);

        return ImmutableQCFailReport.builder()
                .sampleReport(sampleReport)
                .reason(reason)
                .study(study)
                .comments(Optional.ofNullable(comments))
                .isCorrectedReport(correctedReport)
                .signaturePath(reportData.signaturePath())
                .logoRVAPath(reportData.logoRVAPath())
                .logoCompanyPath(reportData.logoCompanyPath())
                .build();
    }
}
