package com.hartwig.hmftools.common.cosmic.fusions;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.EmptyFileException;

import org.junit.Test;

public class COSMICGeneFusionsTest {

    private static final String FUSION_EXAMPLE_FILE = Resources.getResource("cosmic").getPath() + File.separator + "FusionsExample.csv";

    @Test
    public void canReadFromCSV() throws IOException, EmptyFileException {
        COSMICGeneFusionModel model = COSMICGeneFusions.readFromCSV(FUSION_EXAMPLE_FILE);
        assertNotNull(model);
    }
}
