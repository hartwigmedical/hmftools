package com.hartwig.hmftools.patientreporter.qcfail;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;
import java.util.Optional;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumor;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumorFunctions;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.peach.PeachGenotypeFile;
import com.hartwig.hmftools.common.pipeline.PipelineVersionFile;
import com.hartwig.hmftools.common.purple.CheckPurpleQuality;
import com.hartwig.hmftools.common.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContextFile;
import com.hartwig.hmftools.patientreporter.PatientReporterApplication;
import com.hartwig.hmftools.patientreporter.SampleMetadata;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.SampleReportFactory;
import com.hartwig.hmftools.patientreporter.pipeline.PipelineVersion;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class QCFailReporter {

    private static final Logger LOGGER = LogManager.getLogger(QCFailReporter.class);

    @NotNull
    private final QCFailReportData reportData;
    private final String reportDate;

    public QCFailReporter(@NotNull final QCFailReportData reportData, @NotNull final String reportDate) {
        this.reportData = reportData;
        this.reportDate = reportDate;
    }

    @NotNull
    public QCFailReport run(@Nullable QCFailReason reason, @NotNull SampleMetadata sampleMetadata, @NotNull String purplePurityTsv,
            @NotNull String purpleQCFile, @Nullable String comments, boolean correctedReport, @NotNull String expectedPipelineVersion,
            boolean overridePipelineVersion, @Nullable String pipelineVersionFile, boolean requirePipelineVersionFile,
            @NotNull String peachGenotypeTsv) throws IOException {
        assert reason != null;

        String patientId = reportData.limsModel().patientId(sampleMetadata.tumorSampleBarcode());

        PatientPrimaryTumor patientPrimaryTumor =
                PatientPrimaryTumorFunctions.findPrimaryTumorForPatient(reportData.patientPrimaryTumors(), patientId);
        SampleReport sampleReport = SampleReportFactory.fromLimsModel(sampleMetadata, reportData.limsModel(), patientPrimaryTumor);

        if (reason.equals(QCFailReason.SUFFICIENT_TCP_QC_FAILURE) || reason.equals(QCFailReason.INSUFFICIENT_TCP_DEEP_WGS)) {

            String pipelineVersion = null;
            if (requirePipelineVersionFile) {
                assert pipelineVersionFile != null;
                pipelineVersion = PipelineVersionFile.majorDotMinorVersion(pipelineVersionFile);
                PipelineVersion.checkPipelineVersion(pipelineVersion, expectedPipelineVersion, overridePipelineVersion);
            }
        }

        LimsCohortConfig cohort = sampleReport.cohort();

        if (cohort.cohortId().isEmpty()) {
            throw new IllegalStateException("QC fail report not supported for non-cancer study samples: " + sampleMetadata.tumorSampleId());
        }

        String wgsPurityString = null;
        Set<PurpleQCStatus> purpleQc = Sets.newHashSet();
        if (reason.isDeepWGSDataAvailable()) {
            LOGGER.info("Loading PURPLE data from {}", new File(purplePurityTsv).getParent());
            PurityContext purityContext = PurityContextFile.readWithQC(purpleQCFile, purplePurityTsv);

            String formattedPurity = new DecimalFormat("#'%'").format(purityContext.bestFit().purity() * 100);
            boolean hasReliablePurity = CheckPurpleQuality.checkHasReliablePurity(purityContext);

            wgsPurityString = hasReliablePurity ? formattedPurity : Lims.PURITY_NOT_RELIABLE_STRING;
            purpleQc = purityContext.qc().status();
        }

        LOGGER.info("  QC status: {}", purpleQc.toString());

        List<PeachGenotype> peachGenotypesOverrule = Lists.newArrayList();
        if (reason.isDeepWGSDataAvailable() && !purpleQc.contains(PurpleQCStatus.FAIL_CONTAMINATION)) {
            List<PeachGenotype> peachGenotypes = loadPeachData(peachGenotypeTsv);
            peachGenotypesOverrule = sampleReport.reportPharmogenetics() ? peachGenotypes : Lists.newArrayList();
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
                .udiDi(reportData.udiDi())
                .peachGenotypes(peachGenotypesOverrule)
                .purpleQC(purpleQc)
                .reportDate(reportDate)
                .build();
    }

    @NotNull
    private static List<PeachGenotype> loadPeachData(@NotNull String peachGenotypeTsv) throws IOException {
        if (!peachGenotypeTsv.equals(Strings.EMPTY)) {
            LOGGER.info("Loading peach genotypes from {}", new File(peachGenotypeTsv).getParent());
            List<PeachGenotype> peachGenotypes = PeachGenotypeFile.read(peachGenotypeTsv);
            LOGGER.info(" Loaded {} reportable genotypes from {}", peachGenotypes.size(), peachGenotypeTsv);
            return peachGenotypes;
        }else {
            return Lists.newArrayList();
        }
    }
}
