package com.hartwig.hmftools.retentionchecker;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import org.junit.Before;
import org.junit.Test;

import java.net.URL;

import static org.junit.Assert.assertTrue;

public class RetentionCheckAlgoTest {

    private RetentionCheckAlgo algo;

    @Before
    public void setup() {
        algo = new RetentionCheckAlgo(new DummyConverter());
    }

    @Test
    public void algoRunsOnRealTestData() {
        URL originalFastqURL = Resources.getResource("fastq/original");
        String originalFastqPath = originalFastqURL.getPath();

        URL recreatedFastqURL = Resources.getResource("fastq/recreated");
        String recreatedFastqPath = recreatedFastqURL.getPath();

        assertTrue(algo.runAlgo(Lists.newArrayList(originalFastqPath), recreatedFastqPath, 100000));
    }

    @Test
    public void algoRunsOnPartialData() {
        URL originalFastqURL = Resources.getResource("fastq/original");
        String originalFastqPath = originalFastqURL.getPath();

        URL recreatedFastqURL = Resources.getResource("fastq/recreated");
        String recreatedFastqPath = recreatedFastqURL.getPath();

        // KODU: Assume that test files have 25k records.
        assertTrue(algo.runAlgo(Lists.newArrayList(originalFastqPath), recreatedFastqPath, 1000));
    }

    private static class DummyConverter implements FileNameConverter {

        public String apply(String s) {
            return s;
        }
    }
}
