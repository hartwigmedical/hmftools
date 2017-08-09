package com.hartwig.hmftools.patientreportmailer;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import javax.mail.MessagingException;
import javax.mail.Multipart;
import javax.xml.stream.XMLStreamException;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.centra.Centra;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class PatientReportMailerApplication {

    private static final Logger LOGGER = LogManager.getLogger(PatientReportMailerApplication.class);

    private static final String REPORT = "report";
    private static final String DRUP_EMAIL = "drup_email";
    private static final String CENTRA_FILE_PATH = "centra";
    private static final String EMAIL_TEMPLATE = "template";
    private static final String SENDER = "sender";

    public static void main(final String... args)
            throws ParseException, IOException, XMLStreamException, HartwigException, MessagingException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        final String reportPath = cmd.getOptionValue(REPORT);
        final String centraPath = cmd.getOptionValue(CENTRA_FILE_PATH);
        final String templatePath = cmd.getOptionValue(EMAIL_TEMPLATE);
        final String drupEmail = cmd.getOptionValue(DRUP_EMAIL);
        final String sender = cmd.getOptionValue(SENDER);

        if (reportPath == null || centraPath == null || (reportPath.startsWith("DRUP") && drupEmail == null) || sender == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Patient Report Mailer", options);
            System.exit(1);
        }
        final Multipart messageBody = ReportMailer.createMessageBody(templatePath, reportPath);
        final String recipients = getRecipients(reportPath, centraPath, drupEmail);
        final String subject = getSampleFromReportPath(reportPath) + " HMF Report";
        ReportMailer.sendEmail(subject, messageBody, sender, recipients);
        LOGGER.info("Done.");
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(REPORT, true, "Path towards the report");
        options.addOption(DRUP_EMAIL, true, "Drup email (only for DRUP reports)");
        options.addOption(CENTRA_FILE_PATH, true, "Path towards the cpct centra file");
        options.addOption(EMAIL_TEMPLATE, true, "Path towards the email message template");
        options.addOption(SENDER, true, "Sender email address");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    @VisibleForTesting
    static String getRecipients(@NotNull final String reportPath, @NotNull final String centraPath, @NotNull final String drupEmail)
            throws IOException, EmptyFileException {
        final String sample = getSampleFromReportPath(reportPath);
        if (sample.startsWith("DRUP")) {
            final String centraId = getCentraFromSample(sample);
            final Map<String, String> centraRecipients = Centra.readDRUPRecipientsFromCSV(centraPath);
            return centraRecipients.get(centraId).replaceAll(";", ",") + "," + drupEmail;
        } else if (sample.startsWith("CPCT")) {
            final String centraId = getCentraFromSample(sample);
            final Map<String, String> centraRecipients = Centra.readCPCTRecipientsFromCSV(centraPath);
            return centraRecipients.get(centraId).replaceAll(";", ",");
        } else {
            throw new IllegalArgumentException("PatientId was neither CPCT nor DRUP");
        }
    }

    // MIVO: assumes file name is: sample_hmf_report.pdf
    @NotNull
    private static String getSampleFromReportPath(@NotNull final String reportPath) {
        return new File(reportPath).getName().split("_")[0];
    }

    // MIVO: assumes patientId is in DRUP/CPCT format
    @NotNull
    private static String getCentraFromSample(@NotNull final String patientId) {
        return patientId.substring(6, 8);
    }
}
