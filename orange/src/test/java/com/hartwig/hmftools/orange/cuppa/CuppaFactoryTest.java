package com.hartwig.hmftools.orange.cuppa;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.orange.OrangeTestFactory;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class CuppaFactoryTest {

    @Test
    public void canReadFromTestFile() throws IOException {
        List<CuppaEntry> entries = CuppaDataFile.read(OrangeTestFactory.createTestOrangeConfig().cuppaResultCsv());

        CuppaData cuppa = CuppaFactory.build(Strings.EMPTY, entries);

        assertEquals(10, cuppa.simpleDups32To200B());
        assertEquals(8, cuppa.maxComplexSize());
        assertEquals(5, cuppa.LINECount());
        assertEquals(0, cuppa.telomericSGLs());
    }
}