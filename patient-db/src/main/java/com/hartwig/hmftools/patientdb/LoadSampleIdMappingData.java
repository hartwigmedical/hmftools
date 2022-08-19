package com.hartwig.hmftools.patientdb;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.amber.*;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelector;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

public class LoadSampleIdMappingData {

    private static final Logger LOGGER = LogManager.getLogger(LoadSampleIdMappingData.class);

    private static final String SAMPLE = "sample";
    private static final String BARCODE = "barcode";

    public static void main(@NotNull String[] args) throws ParseException, IOException, SQLException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        String sampleId = cmd.getOptionValue(SAMPLE);
        String sampleBarcode = cmd.getOptionValue(BARCODE);

        LOGGER.info("Writing data");
        try (DatabaseAccess dbAccess = databaseAccess(cmd)) {
            dbAccess.writeSampleIdMapping(sampleId, sampleBarcode);
        }
        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Name of tumor sample");
        options.addOption(BARCODE, true, "Barcode of tumor sample");
        addDatabaseCmdLineArgs(options);
        return options;
    }
}
