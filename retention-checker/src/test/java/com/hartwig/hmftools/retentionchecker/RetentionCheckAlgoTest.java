package com.hartwig.hmftools.retentionchecker;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.net.URL;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class RetentionCheckAlgoTest {

    private static final String BASE_FASTQ_PATH = "fastq";

    private final RetentionCheckAlgo algo = new RetentionCheckAlgo(new DummyConverter());

    @Test
    public void algoRunsOnRealTestData() {
        final URL originalFastqURL = resourceURL("original");
        final String originalFastqPath = originalFastqURL.getPath();

        final URL recreatedFastqURL = resourceURL("recreated");
        final String recreatedFastqPath = recreatedFastqURL.getPath();

        assertTrue(algo.runAlgo(Lists.newArrayList(originalFastqPath), recreatedFastqPath, 100000));
    }

    @Test
    public void algoRunsOnPartialData() {
        final URL originalFastqURL = resourceURL("original");
        final String originalFastqPath = originalFastqURL.getPath();

        final URL recreatedFastqURL = resourceURL("recreated");
        final String recreatedFastqPath = recreatedFastqURL.getPath();

        // KODU: Assume that test files have 25k records.
        assertTrue(algo.runAlgo(Lists.newArrayList(originalFastqPath), recreatedFastqPath, 1000));
    }

    @NotNull
    private static URL resourceURL(@NotNull final String fastqPath) {
        return Resources.getResource(BASE_FASTQ_PATH + File.separator + fastqPath);
    }

    private static class DummyConverter implements FileNameConverter {
        public String apply(String s) {
            return s;
        }
    }
}
