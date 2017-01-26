package com.hartwig.hmftools.retentionchecker;

import static org.junit.Assert.assertTrue;

import java.net.URL;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;

import org.junit.Test;

public class RetentionCheckAlgoTest {

    private final RetentionCheckAlgo algo = new RetentionCheckAlgo(new DummyConverter());

    @Test
    public void algoRunsOnRealTestData() {
        final URL originalFastqURL = Resources.getResource("fastq/original");
        final String originalFastqPath = originalFastqURL.getPath();

        final URL recreatedFastqURL = Resources.getResource("fastq/recreated");
        final String recreatedFastqPath = recreatedFastqURL.getPath();

        assertTrue(algo.runAlgo(Lists.newArrayList(originalFastqPath), recreatedFastqPath, 100000));
    }

    @Test
    public void algoRunsOnPartialData() {
        final URL originalFastqURL = Resources.getResource("fastq/original");
        final String originalFastqPath = originalFastqURL.getPath();

        final URL recreatedFastqURL = Resources.getResource("fastq/recreated");
        final String recreatedFastqPath = recreatedFastqURL.getPath();

        // KODU: Assume that test files have 25k records.
        assertTrue(algo.runAlgo(Lists.newArrayList(originalFastqPath), recreatedFastqPath, 1000));
    }

    private static class DummyConverter implements FileNameConverter {
        public String apply(String s) {
            return s;
        }
    }
}
