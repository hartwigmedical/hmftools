package com.hartwig.hmftools.common.copynumber.freec;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.HartwigException;

import org.junit.Test;

public class FreecGCContentFactoryTest {

    private static final String BASE_PATH = Resources.getResource("copynumber").getPath() + File.separator + "gccontent";

    @Test
    public void canLoadNormalFile() throws IOException, HartwigException {
        final List<FreecGCContent> gcContent = FreecGCContentFactory.loadGCContent(BASE_PATH);
        assertEquals(100, gcContent.size());
    }
}
