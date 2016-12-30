package com.hartwig.healthchecker.runners;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.google.common.io.Resources;
import com.hartwig.healthchecker.runners.checks.HealthCheck;
import com.hartwig.healthchecker.runners.checks.MetadataCheck;
import com.hartwig.healthchecker.exception.HealthChecksException;
import com.hartwig.healthchecker.io.dir.RunContext;
import com.hartwig.healthchecker.io.dir.TestRunContextFactory;
import com.hartwig.healthchecker.result.BaseResult;
import com.hartwig.healthchecker.result.PatientResult;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class MetadataCheckerTest {

    private static final String RUN_DIRECTORY = Resources.getResource("checks/metadata").getPath();

    private static final String EXPECTED_VERSION = "v1.7";
    private static final String EXPECTED_DATE = "2016-07-09";
    private static final int EXPECTED_NUM_CHECKS = 4;

    private final MetadataChecker checker = new MetadataChecker();

    @Test
    public void canReadCorrectMetaData() throws IOException, HealthChecksException {
        final RunContext runContext = TestRunContextFactory.forTest(RUN_DIRECTORY);
        final BaseResult result = checker.tryRun(runContext);

        final PatientResult patientResult = (PatientResult) result;
        assertEquals(EXPECTED_NUM_CHECKS, patientResult.getRefSampleChecks().size());
        assertEquals(EXPECTED_NUM_CHECKS, patientResult.getTumorSampleChecks().size());
        assertCheck(patientResult.getRefSampleChecks(), MetadataCheck.RUN_DATE, EXPECTED_DATE);
        assertCheck(patientResult.getRefSampleChecks(), MetadataCheck.PIPELINE_VERSION, EXPECTED_VERSION);
        assertCheck(patientResult.getTumorSampleChecks(), MetadataCheck.RUN_DATE, EXPECTED_DATE);
        assertCheck(patientResult.getTumorSampleChecks(), MetadataCheck.PIPELINE_VERSION, EXPECTED_VERSION);
    }

    @Test
    public void errorYieldsCorrectOutput() {
        final RunContext runContext = TestRunContextFactory.forTest(RUN_DIRECTORY);
        final PatientResult result = (PatientResult) checker.errorRun(runContext);
        assertEquals(EXPECTED_NUM_CHECKS, result.getRefSampleChecks().size());
        assertEquals(EXPECTED_NUM_CHECKS, result.getTumorSampleChecks().size());
    }

    private static void assertCheck(@NotNull final List<HealthCheck> checks, @NotNull final MetadataCheck checkName,
            @NotNull final String expectedValue) {
        final Optional<HealthCheck> optCheck = checks.stream().filter(
                check -> check.getCheckName().equals(checkName.toString())).findFirst();
        assert optCheck.isPresent();

        assertEquals(expectedValue, optCheck.get().getValue());
    }
}