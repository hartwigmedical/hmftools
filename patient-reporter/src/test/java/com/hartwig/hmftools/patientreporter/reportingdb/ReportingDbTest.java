package com.hartwig.hmftools.patientreporter.reportingdb;

import static java.lang.String.join;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Map;

import com.google.common.io.Resources;
import com.google.gson.GsonBuilder;
import com.google.gson.reflect.TypeToken;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortTestFactory;
import com.hartwig.hmftools.common.lims.reportingdb.ReportingDatabase;
import com.hartwig.hmftools.common.lims.reportingdb.ReportingEntry;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.patientreporter.ExampleAnalysisConfig;
import com.hartwig.hmftools.patientreporter.ExampleAnalysisTestFactory;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;

import org.junit.Test;

public class ReportingDbTest {

    @Test
    public void shouldWriteApiUpdateJson() throws IOException {
        File reportingDbTsv = File.createTempFile("reporting-db-test", ".tsv");

        BufferedWriter writer = new BufferedWriter(new FileWriter(reportingDbTsv, true));
        writer.write("tumorBarcode\tsampleId\tcohort\treportDate\treportType\tpurity\thasReliableQuality\thasReliablePurity\n");
        writer.close();
        LimsCohortConfig cohortConfig = LimsCohortTestFactory.createCPCTCohortConfig();

        ExampleAnalysisConfig config =
                new ExampleAnalysisConfig.Builder().sampleId("CPCT01_SUCCESS").limsCohortConfig(cohortConfig).build();

        ReportingDb reportingDb = new ReportingDb();

        File expectedOutput = new File("/tmp/CPCT01_SUCCESS_FR12345678_dna_analysis_report_api-update.json");
        Files.deleteIfExists(expectedOutput.toPath());
        assertFalse(expectedOutput.exists());
        reportingDb.appendAnalysedReport(ExampleAnalysisTestFactory.createAnalysisWithAllTablesFilledIn(config, PurpleQCStatus.PASS),
                "/tmp");
        assertTrue(expectedOutput.exists());

        Map<String, Object> output = new GsonBuilder().serializeNulls()
                .serializeSpecialFloatingPointValues()
                .create()
                .fromJson(join("\n", Files.readAllLines(expectedOutput.toPath())), new TypeToken<Map<String, Object>>() {
                }.getType());
        assertEquals(output.get("has_reliable_quality"), true);
        assertEquals(output.get("has_reliable_purity"), true);
        assertEquals(output.get("purity"), 1.0);
        assertEquals(output.get("cohort"), "CPCT");
        assertEquals(output.get("report_type"), "dna_analysis_report");
        assertEquals(output.get("barcode"), "FR12345678");
        assertEquals(output.get("report_date"), LocalDateTime.now().format(DateTimeFormatter.ofPattern("dd-MMM-y")));
    }
}