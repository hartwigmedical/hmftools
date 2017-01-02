package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.hartwig.hmftools.common.exception.HealthChecksException;
import com.hartwig.hmftools.healthchecker.context.RunContext;
import com.hartwig.hmftools.healthchecker.context.TestRunContextFactory;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.SingleValueResult;
import com.hartwig.hmftools.healthchecker.runners.checks.HealthCheck;
import com.hartwig.hmftools.healthchecker.runners.checks.SlicedCheck;

import org.junit.Assert;
import org.junit.Test;

public class SlicedCheckerTest {

    private static final String RUN_DIRECTORY = RunnerTestFunctions.getRunnerResourcePath("sliced");
    private static final String REF_SAMPLE = "sample1";
    private static final String TUMOR_SAMPLE = "sample2";

    private static final int EXPECTED_SLICED_COUNT = 4;

    private final SlicedChecker checker = new SlicedChecker();

    @Test
    public void canAnalyseTypicalSlicedVCF() throws IOException, HealthChecksException {
        final RunContext runContext = TestRunContextFactory.forTest(RUN_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);

        final BaseResult result = checker.tryRun(runContext);
        Assert.assertEquals(CheckType.SLICED, result.getCheckType());

        final HealthCheck check = ((SingleValueResult) result).getCheck();
        Assert.assertEquals(SlicedCheck.SLICED_NUMBER_OF_VARIANTS.toString(), check.getCheckName());
        assertEquals(TUMOR_SAMPLE, check.getSampleId());
        assertEquals(Integer.toString(EXPECTED_SLICED_COUNT), check.getValue());
    }

    @Test
    public void errorYieldsCorrectOutput() {
        final RunContext runContext = TestRunContextFactory.forTest(RUN_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);
        final SingleValueResult result = (SingleValueResult) checker.errorRun(runContext);
        assertNotNull(result.getCheck());
    }

    @Test(expected = IOException.class)
    public void readingNonExistingFileYieldsIOException() throws IOException, HealthChecksException {
        final RunContext runContext = TestRunContextFactory.forTest("DoesNotExist", REF_SAMPLE, TUMOR_SAMPLE);
        checker.tryRun(runContext);
    }
}
