package com.hartwig.hmftools.sullivan;

import com.google.common.io.Resources;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.net.URL;
import java.io.File;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class SullivanAlgoTest {

    private SullivanAlgo algo;

    @Before
    public void setup() {
        algo = new SullivanAlgo(new DummyConverter());
    }
    @Test
    public void sullivanAlgoRunsOnRealTestData() {
        URL originalFastqURL = Resources.getResource("fastq/original.fastq");
        String originalFastqPath = originalFastqURL.getPath();

        URL recreatedFastqURL = Resources.getResource("fastq/recreated.fastq");
        String recreatedFastqPath = recreatedFastqURL.getPath();

        assertTrue(algo.runSullivanAlgo(originalFastqPath, recreatedFastqPath, null, false, 100000));
    }

    @Test
    public void sullivanAlgoRunsOnPartialData() {
        URL originalFastqURL = Resources.getResource("fastq/original.fastq");
        String originalFastqPath = originalFastqURL.getPath();

        URL recreatedFastqURL = Resources.getResource("fastq/recreated.fastq");
        String recreatedFastqPath = recreatedFastqURL.getPath();

        // KODU: Assume that test files have 25k records.
        assertTrue(algo.runSullivanAlgo(originalFastqPath, recreatedFastqPath, null, false, 1000));
    }

    // KODU: Should mock this.
    private static class DummyConverter implements FileNameConverter {

        public String apply(String s) {
            return s;
        }
    }

}
