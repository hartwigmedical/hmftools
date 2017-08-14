package com.hartwig.hmftools.patientreportmailer;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.EmptyFileException;

import org.junit.Test;

public class PatientReportMailerTest {

    private static final String CENTRA_FILE = Resources.getResource("centra.csv").getPath();
    private static final String CPCT_REPORT = Resources.getResource("CPCT02010001T_hmf_report.pdf").getPath();
    private static final String DRUP_REPORT = Resources.getResource("DRUP01010001TII_hmf_report.pdf").getPath();
    private static final String OTHER_REPORT = Resources.getResource("TEST_hmf_report.pdf").getPath();
    private static final String TEMPLATE_FILE = Resources.getResource("template.txt").getPath();
    private static final String MEB_DATE = "03-03-2017";
    private static final String MEB_DEADLINE = "02-02-2017";
    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("dd-MM-yyyy");

    @Test
    public void canReadCPCTRecipients() throws IOException, EmptyFileException {
        final String recipients = PatientReportMailerApplication.getRecipients(CPCT_REPORT, CENTRA_FILE, "");
        assertEquals("my@email.com, my2@email.com", recipients);
    }

    @Test
    public void canReadDRUPRecipients() throws IOException, EmptyFileException {
        final String recipients = PatientReportMailerApplication.getRecipients(DRUP_REPORT, CENTRA_FILE, "drup@email.com");
        assertEquals("my3@email.com,drup@email.com", recipients);
    }

    @Test(expected = IllegalArgumentException.class)
    public void throwsOnNonCPCTorDRUP() throws IOException, EmptyFileException {
        PatientReportMailerApplication.getRecipients(OTHER_REPORT, CENTRA_FILE, "drup@email.com");
    }

    @Test
    public void correctlyFillsTemplate() throws IOException, EmptyFileException {
        final LocalDate mebDate = LocalDate.parse(MEB_DATE, DATE_FORMATTER);
        final LocalDate mebDeadline = LocalDate.parse(MEB_DEADLINE, DATE_FORMATTER);
        final String messageBody = ReportMailer.fillTemplate(TEMPLATE_FILE, mebDate, mebDeadline);
        final String[] messageLines = messageBody.split("\n");
        assertEquals("meb.date: " + MEB_DATE, messageLines[0]);
        assertEquals("meb.deadline: " + MEB_DEADLINE, messageLines[1]);
    }

    @Test
    public void determinesNextMebDates() {
        final LocalDate mebSeed = LocalDate.parse("15-08-2017", DATE_FORMATTER);
        final LocalDate expectedDate = LocalDate.parse("10-08-2017", DATE_FORMATTER);
        assertEquals(expectedDate, ReportMailer.determineMebDeadline(mebSeed, expectedDate));
        assertEquals(LocalDate.parse("24-08-2017", DATE_FORMATTER), ReportMailer.determineMebDeadline(mebSeed, expectedDate.plusDays(1)));
    }
}
