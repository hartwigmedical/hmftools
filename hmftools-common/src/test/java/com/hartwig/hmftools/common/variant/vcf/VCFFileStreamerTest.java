package com.hartwig.hmftools.common.variant.vcf;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.io.BufferedReader;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.GermlineVariant;

import org.junit.Test;

public class VCFFileStreamerTest {

    private static final String RESOURCE_DIR = Resources.getResource("variants").getPath();

    @Test
    public void canStreamGermlineVCFAsExpected() throws IOException, HartwigException {
        final BufferedReader reader = VCFFileStreamer.getVCFReader(RESOURCE_DIR, "germline.vcf");

        final GermlineVariant firstVariant = VCFFileStreamer.nextVariant(reader);
        assertNotNull(firstVariant);
        assertEquals("0/1:12,17:29:99:500,0,440", firstVariant.refData());

        // KODU: There are 2 more variants in the test file.
        assertNotNull(VCFFileStreamer.nextVariant(reader));
        assertNotNull(VCFFileStreamer.nextVariant(reader));

        assertNull(VCFFileStreamer.nextVariant(reader));
        // KODU: Once null -> remain null
        assertNull(VCFFileStreamer.nextVariant(reader));
    }
}