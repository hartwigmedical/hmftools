package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.sigs.SignatureAllocation;
import com.hartwig.hmftools.common.sigs.SignatureAllocationFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadSignatures {
    private static final Logger LOGGER = LogManager.getLogger(LoadSignatures.class);

    private static final String SAMPLE = "sample";
    private static final String ISOLATION_BARCODE = "isolation_barcode";
    private static final String SAMPLE_DIR = "sample_dir";

    public static void main(@NotNull String[] args) throws ParseException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);
        DatabaseAccess dbAccess = createDatabaseAccess(cmd);

        if (dbAccess == null) {

            LOGGER.error("Failed to create DB connection");
            System.exit(1);
        }

        String sampleId = cmd.getOptionValue(SAMPLE);
        String isolationBarcode = cmd.getOptionValue(ISOLATION_BARCODE);
        String sampleDir = cmd.getOptionValue(SAMPLE_DIR);

        loadSignatureData(dbAccess, sampleId, isolationBarcode, sampleDir);

        LOGGER.info("signature allocation loading complete");
    }

    private static void loadSignatureData(final DatabaseAccess dbAccess, final String sampleId, @NotNull String isolationBarcode,
            final String sampleDir) {
        try {
            final List<SignatureAllocation> sigAllocations =
                    SignatureAllocationFile.read(SignatureAllocationFile.generateFilename(sampleDir, sampleId));

            if (!sigAllocations.isEmpty()) {
                LOGGER.info("sample({}) writing {} allocations to database", sampleId, sigAllocations.size());
                dbAccess.writeSignatures(sampleId, isolationBarcode, sigAllocations);
            } else {
                LOGGER.info("sample({}) has not signature allocations", sampleId);
            }
        } catch (IOException e) {
            LOGGER.error("failed to load sample({}) allocations: {}", sampleId, e.toString());
        }
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        addDatabaseCmdLineArgs(options);
        options.addOption(SAMPLE, true, "Name of the tumor sample");
        options.addOption(ISOLATION_BARCODE, true, "Name of the isolation barcode");
        options.addOption(SAMPLE_DIR, true, "Directory to read signature data from");

        return options;
    }
}
