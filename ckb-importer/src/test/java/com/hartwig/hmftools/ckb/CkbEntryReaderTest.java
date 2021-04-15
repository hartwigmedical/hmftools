package com.hartwig.hmftools.ckb;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class CkbEntryReaderTest {

    private static final String CKB_DIR = Resources.getResource("ckb").getPath();

    @Test
    public void canReadFromTestDir() throws IOException {
        assertNotNull(CkbEntryReader.read(CKB_DIR));
    }
}