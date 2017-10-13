package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.context.TestRunContextFactory;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;

import org.jetbrains.annotations.NotNull;
import org.junit.Assert;
import org.junit.Test;

public class PurpleCheckerTest {
    private static final String BASE_DIRECTORY = Resources.getResource("purple").getPath();
    private static final String CORRECT_RUN = BASE_DIRECTORY + File.separator + "run";
    private static final String REF_SAMPLE = "refSample";
    private static final String TUMOR_SAMPLE = "tumorSample";

    private final PurpleChecker checker = new PurpleChecker();

    @Test
    public void extractDataFromPurpleWorksForSomatic() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(CORRECT_RUN, REF_SAMPLE, TUMOR_SAMPLE);
        final BaseResult result = checker.tryRun(runContext);

        Assert.assertEquals(CheckType.PURPLE, result.getCheckType());
        final List<HealthCheck> checks = ((MultiValueResult) result).getChecks();

        assertCheck(checks, PurpleCheck.PURPLE_GENDER.toString(), "MALE:MALE");
        assertCheck(checks, PurpleCheck.PURPLE_SEGMENT_SCORE.toString(), "88");
    }


    private static void assertCheck(@NotNull final List<HealthCheck> checks, @NotNull final String checkName, final String expectedValue) {
        final Optional<HealthCheck> report = checks.stream().filter(data -> data.getCheckName().equals(checkName)).findFirst();

        assert report.isPresent();
        final String check = report.get().getValue();
        assertEquals(expectedValue, check);
    }

}
