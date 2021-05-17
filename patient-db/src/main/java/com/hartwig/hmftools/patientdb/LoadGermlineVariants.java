package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.sql.SQLException;

import com.hartwig.hmftools.patientdb.dao.BufferedWriter;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class LoadGermlineVariants {

    private static final Logger LOGGER = LogManager.getLogger(LoadGermlineVariants.class);

    private static final String SAMPLE = "sample";
    private static final String REFERENCE = "reference";
    private static final String RNA = "rna";

    private static final String GERMLINE_VCF = "germline_vcf";

    public static void main(@NotNull String[] args) throws ParseException, IOException, SQLException {
        Options options = createBasicOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        String vcfFileLocation = cmd.getOptionValue(GERMLINE_VCF);
        String sample = cmd.getOptionValue(SAMPLE);
        String referenceSample = cmd.getOptionValue(REFERENCE, Strings.EMPTY);
        String rnaSample = cmd.getOptionValue(RNA, Strings.EMPTY);
        DatabaseAccess dbAccess = databaseAccess(cmd);

        LOGGER.info("Removing old data of sample {}", sample);
        try (AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFileLocation, new VCFCodec(), false);
                BufferedWriter<VariantContext> dbWriter = dbAccess.germlineVariantWriter(sample, referenceSample, rnaSample)) {
            LOGGER.info("Streaming data from {} to db", vcfFileLocation);
            dbWriter.initialise();

            for (VariantContext context : reader.iterator()) {
                dbWriter.accept(context);
            }
        }

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createBasicOptions() {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Name of the tumor sample.");
        options.addOption(REFERENCE, true, "Optional name of the reference sample. ");
        options.addOption(RNA, true, "Optional name of the rna sample.");

        options.addOption(GERMLINE_VCF, true, "Path to the germline SNV/indel vcf file.");
        addDatabaseCmdLineArgs(options);

        return options;
    }
}
