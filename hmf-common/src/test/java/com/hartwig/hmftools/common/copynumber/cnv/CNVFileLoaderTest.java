package com.hartwig.hmftools.common.copynumber.cnv;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.exception.HartwigException;

import org.junit.Test;

public class CNVFileLoaderTest {
    private static final String BASE_PATH = Resources.getResource("copynumber").getPath();
    private static final String SAMPLE = "sample";

    @Test
    public void canLoadCNVFile() throws IOException, HartwigException {
        final List<CopyNumber> copyNumbers = CNVFileLoader.loadCNV(BASE_PATH, SAMPLE);
        assertEquals(2, copyNumbers.size());
    }
}