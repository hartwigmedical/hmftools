package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.healthchecker.context.RunContext;
import com.hartwig.hmftools.healthchecker.context.TestRunContextFactory;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;
import com.hartwig.hmftools.healthchecker.result.PatientResult;
import com.hartwig.hmftools.healthchecker.runners.checks.HealthCheck;
import com.hartwig.hmftools.healthchecker.runners.checks.MetadataCheck;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class MetadataCheckerTest {

    private static final String METADATA_BASE_PATH = RunnerTestFunctions.getRunnerResourcePath("metadata");

    private static final String RUN_DIRECTORY_SOMATIC = METADATA_BASE_PATH + File.separator + "run1";
    private static final String EXPECTED_VERSION_SOMATIC = "v1.7";
    private static final String EXPECTED_DATE_SOMATIC = "2016-07-09";

    private static final String RUN_DIRECTORY_SINGLE_SAMPLE = METADATA_BASE_PATH + File.separator + "run2";
    private static final String EXPECTED_VERSION_SINGLE_SAMPLE = "v1.11-rc.2";
    private static final String EXPECTED_DATE_SINGLE_SAMPLE = "2017-01-01";

    private static final int EXPECTED_NUM_CHECKS = 4;
    private static final String REF_SAMPLE = "sample1";
    private static final String TUMOR_SAMPLE = "sample2";

    private final MetadataChecker checker = new MetadataChecker();

    @Test
    public void canReadCorrectMetaDataSomatic() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY_SOMATIC, REF_SAMPLE,
                TUMOR_SAMPLE);
        final BaseResult result = checker.tryRun(runContext);

        final PatientResult patientResult = (PatientResult) result;
        assertEquals(EXPECTED_NUM_CHECKS, patientResult.getRefSampleChecks().size());
        assertEquals(EXPECTED_NUM_CHECKS, patientResult.getTumorSampleChecks().size());
        assertCheck(patientResult.getRefSampleChecks(), MetadataCheck.RUN_DATE, EXPECTED_DATE_SOMATIC);
        assertCheck(patientResult.getRefSampleChecks(), MetadataCheck.PIPELINE_VERSION, EXPECTED_VERSION_SOMATIC);
        assertCheck(patientResult.getTumorSampleChecks(), MetadataCheck.RUN_DATE, EXPECTED_DATE_SOMATIC);
        assertCheck(patientResult.getTumorSampleChecks(), MetadataCheck.PIPELINE_VERSION, EXPECTED_VERSION_SOMATIC);
    }

    @Test
    public void canReadCorrectMetaDataSingleSample() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSingleSampleTest(RUN_DIRECTORY_SINGLE_SAMPLE,
                REF_SAMPLE);
        final MultiValueResult result = (MultiValueResult) checker.tryRun(runContext);

        assertEquals(EXPECTED_NUM_CHECKS, result.getChecks().size());
        assertCheck(result.getChecks(), MetadataCheck.RUN_DATE, EXPECTED_DATE_SINGLE_SAMPLE);
        assertCheck(result.getChecks(), MetadataCheck.PIPELINE_VERSION, EXPECTED_VERSION_SINGLE_SAMPLE);
    }

    @Test
    public void errorYieldsCorrectOutputForSomatic() {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY_SOMATIC, REF_SAMPLE,
                TUMOR_SAMPLE);
        final PatientResult result = (PatientResult) checker.errorRun(runContext);
        assertEquals(EXPECTED_NUM_CHECKS, result.getRefSampleChecks().size());
        assertEquals(EXPECTED_NUM_CHECKS, result.getTumorSampleChecks().size());
    }

    @Test
    public void errorYieldsCorrectOutputForSingleSample() {
        final RunContext runContext = TestRunContextFactory.forSingleSampleTest(RUN_DIRECTORY_SINGLE_SAMPLE,
                REF_SAMPLE);
        final MultiValueResult result = (MultiValueResult) checker.errorRun(runContext);
        assertEquals(EXPECTED_NUM_CHECKS, result.getChecks().size());
    }

    private static void assertCheck(@NotNull final List<HealthCheck> checks, @NotNull final MetadataCheck checkName,
            @NotNull final String expectedValue) {
        final Optional<HealthCheck> optCheck = checks.stream().filter(
                check -> check.getCheckName().equals(checkName.toString())).findFirst();
        assert optCheck.isPresent();

        assertEquals(expectedValue, optCheck.get().getValue());
    }
}