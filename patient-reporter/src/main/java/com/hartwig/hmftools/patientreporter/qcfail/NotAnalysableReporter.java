package com.hartwig.hmftools.patientreporter.qcfail;

import java.io.File;
import java.io.IOException;
import java.util.Optional;

import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocationFunctions;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.patientreporter.BaseReportData;
import com.hartwig.hmftools.patientreporter.ImmutableNotAnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.ImmutableSampleReport;
import com.hartwig.hmftools.patientreporter.NotAnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.SampleReport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class NotAnalysableReporter {

    private static final Logger LOGGER = LogManager.getLogger(NotAnalysableReporter.class);
    private static final String PURPLE_DIRECTORY = "purple";

    @NotNull
    abstract BaseReportData baseReportData();

    public NotAnalysedPatientReport run(@NotNull final String sample, @NotNull final NotAnalysableReason reason,
            @Nullable final String comments) throws IOException { // @Nullable final String runDir, final boolean doShallowPurity,
        final Lims lims = baseReportData().limsModel();

        final NotAnalysableStudy study = NotAnalysableStudy.fromSample(sample);

       // final RunContext run = ProductionRunContextFactory.fromRunDirectory(runDir);
       // final double shallowPurity = doShallowPurity ? extractPurity(run) : null;

        assert study != null;

        final PatientTumorLocation patientTumorLocation =
                PatientTumorLocationFunctions.findPatientTumorLocationForSample(baseReportData().patientTumorLocations(), sample);

        final SampleReport sampleReport = ImmutableSampleReport.of(sample,
                patientTumorLocation,
               // lims.labelSample(sample).equalsIgnoreCase("core") ? Double.toString(shallowPurity) : lims.tumorPercentageForSample(sample),
                lims.tumorPercentageForSample(sample),
                lims.arrivalDateForSample(sample),
                null,
                lims.labProceduresForSample(sample),
                lims.labelSample(sample).equalsIgnoreCase("core")
                        ? baseReportData().centerModel()
                        .getCoreRecipients(lims.projectNameDVO(sample))
                        : baseReportData().centerModel().getAddresseeStringForSample(sample),
                lims.labelSample(sample),
                lims.projectNameDVO(sample),
                lims.labelSample(sample).equalsIgnoreCase("core") ? lims.contactEmail(sample) : "",
                lims.labelSample(sample).equalsIgnoreCase("core") ? lims.contactName(sample) : "",
                lims.patientNumber(sample));

        return ImmutableNotAnalysedPatientReport.of(sampleReport,
                reason,
                study,
                Optional.ofNullable(comments),
                baseReportData().signaturePath(),
                baseReportData().logoRVAPath());
    }

    private double extractPurity(@NotNull RunContext runDir) throws IOException {
        LOGGER.info("Loading purple data for sample " + runDir.tumorSample());
        final String cnvBasePath = runDir.runDirectory() + File.separator + PURPLE_DIRECTORY;

        final PurityContext purityContext = FittedPurityFile.read(cnvBasePath, runDir.tumorSample());

        return purityContext.bestFit().purity();
    }
}
