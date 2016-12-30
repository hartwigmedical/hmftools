package com.hartwig.healthchecker.report;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Optional;

import com.hartwig.healthchecker.runners.CheckType;
import com.hartwig.healthchecker.exception.GenerateReportException;
import com.hartwig.healthchecker.exception.HealthChecksException;
import com.hartwig.healthchecker.io.dir.RunContext;
import com.hartwig.healthchecker.io.dir.TestRunContextFactory;
import com.hartwig.healthchecker.result.BaseResult;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import mockit.Mocked;
import mockit.NonStrictExpectations;

public class ReportTest {

    private static final String SOME_VERSION = "v1.7";
    private static final String SOME_DATE = "2016-07-09";

    private static final RunContext MOCK_RUN_CONTEXT = TestRunContextFactory.forTest("path/to/run");
    private static final String MOCK_OUTPUT_DIR = "path/to/output";

    @Test
    public void generateStdOutReport() throws IOException, HealthChecksException {
        final Report report = StandardOutputReport.getInstance();

        final BaseResult baseConfig1 = new TestResult(CheckType.MAPPING);
        report.addResult(baseConfig1);

        final BaseResult baseConfig2 = new TestResult(CheckType.PRESTATS);
        report.addResult(baseConfig2);

        final Optional<String> jsonOptional = report.generateReport(MOCK_RUN_CONTEXT, MOCK_OUTPUT_DIR);
        assertNotNull(jsonOptional);
        assertTrue(jsonOptional.isPresent());
        final String json = jsonOptional.get();
        assertFalse(json.contains(SOME_DATE));
        assertFalse(json.contains(SOME_VERSION));
    }

    @Test
    public void generateReportMetadataIOException() throws IOException, HealthChecksException {
        final Report report = StandardOutputReport.getInstance();

        final BaseResult baseConfig1 = new TestResult(CheckType.MAPPING);
        report.addResult(baseConfig1);

        final BaseResult baseConfig2 = new TestResult(CheckType.PRESTATS);
        report.addResult(baseConfig2);

        final Optional<String> jsonOptional = report.generateReport(MOCK_RUN_CONTEXT, MOCK_OUTPUT_DIR);
        assertNotNull(jsonOptional);
        assertTrue(jsonOptional.isPresent());
        final String json = jsonOptional.get();
        assertFalse(json.contains(SOME_DATE));
        assertFalse(json.contains(SOME_VERSION));
    }

    @Test
    public void generateReportMetadataHealthCheckException() throws IOException, HealthChecksException {
        final Report report = StandardOutputReport.getInstance();

        final BaseResult baseConfig1 = new TestResult(CheckType.MAPPING);
        report.addResult(baseConfig1);

        final BaseResult baseConfig2 = new TestResult(CheckType.PRESTATS);
        report.addResult(baseConfig2);

        final Optional<String> jsonOptional = report.generateReport(MOCK_RUN_CONTEXT, MOCK_OUTPUT_DIR);
        assertNotNull(jsonOptional);
        assertTrue(jsonOptional.isPresent());
        final String json = jsonOptional.get();
        assertFalse(json.contains(SOME_DATE));
        assertFalse(json.contains(SOME_VERSION));
    }

    @Test(expected = GenerateReportException.class)
    public void generateReportException(@Mocked final FileWriter fileWriter)
            throws IOException, HealthChecksException {

        new NonStrictExpectations() {
            {
                new FileWriter(new File(anyString));
                result = fileWriter;
                times = 1;

                fileWriter.write(anyString);
                result = new IOException("");
            }
        };
        final Report report = JsonReport.getInstance();

        final BaseResult baseConfig1 = new TestResult(CheckType.MAPPING);
        report.addResult(baseConfig1);

        final Optional<String> result = report.generateReport(MOCK_RUN_CONTEXT, MOCK_OUTPUT_DIR);

        assertNotNull(result);
        assertFalse(result.isPresent());
    }

    private static class TestResult implements BaseResult {
        @NotNull
        private final CheckType checkType;

        TestResult(@NotNull final CheckType checkType) {
            this.checkType = checkType;
        }

        @NotNull
        @Override
        public CheckType getCheckType() {
            return checkType;
        }
    }
}
