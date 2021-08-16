package com.hartwig.hmftools.orange.cuppa;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.orange.OrangeTestFactory;

import org.junit.Test;

public class CuppaDataFileTest {

    @Test
    public void canLoadTestFile() throws IOException {
        List<CuppaEntry> entries = CuppaDataFile.read(OrangeTestFactory.createTestOrangeConfig().cuppaResultCsv());
        assertEquals(1044, entries.size());
    }
}