package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class LimsWideFileTest {

    private static final String FILE = Resources.getResource("lims/contact_info.tsv").getPath();
    private static final String STUDY_NAME = "AAA";
    private static final String REPORT_RECEIVER_NAME = "contact";
    private static final String REPORT_RECEIVER_EMIAL = "contact.nl";

    @Test
    public void canReadLimsWide() throws IOException {
        LimsWide limsWide = LimsWideFile.read(FILE);

        assertEquals(STUDY_NAME, limsWide.studyName());
        assertEquals(REPORT_RECEIVER_NAME, limsWide.reportReceiverName());
        assertEquals(REPORT_RECEIVER_EMIAL, limsWide.reportReceiverEmail());
    }
}