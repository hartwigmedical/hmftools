package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.context.TestRunContextFactory;
import com.hartwig.hmftools.common.io.exception.MalformedFileException;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;

import org.jetbrains.annotations.NotNull;
import org.junit.Assert;
import org.junit.Test;

public class PurpleCheckerTest {
    private static final String BASE_DIRECTORY = Resources.getResource("").getPath();
    private static final String REF_SAMPLE = "refSample";
    private static final String TUMOR_SAMPLE = "tumorSample";

    private final PurpleChecker checker = new PurpleChecker();

    @Test
    public void extractDataFromPurpleWorksForSomatic() throws IOException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(BASE_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);
        final BaseResult result = checker.run(runContext);

        Assert.assertEquals(CheckType.PURPLE, result.checkType());
        final List<HealthCheck> checks = ((MultiValueResult) result).checks();

        assertCheck(checks, PurpleCheck.AMBER_GENDER.toString(), "MALE");
        assertCheck(checks, PurpleCheck.COBALT_GENDER.toString(), "FEMALE");
        assertCheck(checks, PurpleCheck.PURPLE_SEGMENT_SCORE.toString(), "199");
        assertCheck(checks, PurpleCheck.PURPLE_DELETED_GENES_SCORE.toString(), "120");
    }

    @Test(expected = MalformedFileException.class)
    public void testMalformed() throws IOException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(BASE_DIRECTORY, REF_SAMPLE, "malformed");
        checker.run(runContext);
    }

    @Test(expected = IOException.class)
    public void testMissing() throws IOException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(BASE_DIRECTORY, REF_SAMPLE, "missing");
        checker.run(runContext);
    }

    private static void assertCheck(@NotNull final List<HealthCheck> checks, @NotNull final String checkName, final String expectedValue) {
        final Optional<HealthCheck> report = checks.stream().filter(data -> data.getCheckName().equals(checkName)).findFirst();

        assert report.isPresent();
        final String check = report.get().getValue();
        assertEquals(expectedValue, check);
    }
}
