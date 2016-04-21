package com.hartwig.hmftools.sullivan;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import org.junit.Before;
import org.junit.Test;

import java.net.URL;

import static org.junit.Assert.assertTrue;

public class SullivanAlgoTest {

    private SullivanAlgo algo;

    @Before
    public void setup() {
        algo = new SullivanAlgo(new DummyConverter());
    }

    @Test
    public void sullivanAlgoRunsOnRealTestData() {
        URL originalFastqURL = Resources.getResource("fastq/original");
        String originalFastqPath = originalFastqURL.getPath();

        URL recreatedFastqURL = Resources.getResource("fastq/recreated");
        String recreatedFastqPath = recreatedFastqURL.getPath();

        assertTrue(algo.runSullivanAlgo(Lists.newArrayList(originalFastqPath), recreatedFastqPath, 100000));
    }

    @Test
    public void sullivanAlgoRunsOnPartialData() {
        URL originalFastqURL = Resources.getResource("fastq/original");
        String originalFastqPath = originalFastqURL.getPath();

        URL recreatedFastqURL = Resources.getResource("fastq/recreated");
        String recreatedFastqPath = recreatedFastqURL.getPath();

        // KODU: Assume that test files have 25k records.
        assertTrue(algo.runSullivanAlgo(Lists.newArrayList(originalFastqPath), recreatedFastqPath, 1000));
    }

    private static class DummyConverter implements FileNameConverter {

        public String apply(String s) {
            return s;
        }
    }
}
