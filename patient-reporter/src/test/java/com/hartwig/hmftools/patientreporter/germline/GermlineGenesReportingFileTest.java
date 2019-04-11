package com.hartwig.hmftools.patientreporter.germline;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Set;
import java.util.stream.Collectors;

import org.junit.Test;

public class GermlineGenesReportingFileTest {

    @Test
    public void doReadGermlineFilesCorrectly() {
        final InputStream inputStreamGermlineGenesReported =
                GermlineGenesReportingFile.class.getResourceAsStream("/germline/GermlineGenesReported.tsv");
        final InputStream inputStreamGermlineNotify =
                GermlineGenesReportingFile.class.getResourceAsStream("/germline/GermlineVariantNotify.tsv");

        Set<String> germlineGenesReported =
                new BufferedReader(new InputStreamReader(inputStreamGermlineGenesReported)).lines().collect(Collectors.toSet());
        Set<String> germlineVariantNotify =
                new BufferedReader(new InputStreamReader(inputStreamGermlineNotify)).lines().collect(Collectors.toSet());

        assertEquals(germlineGenesReported, GermlineGenesReportingFile.readingLines().germlineGenes());
        assertEquals(germlineVariantNotify, GermlineGenesReportingFile.readingLines().germlineGenesNotify());
    }
}