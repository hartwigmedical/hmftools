package com.hartwig.hmftools.patientdb;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.bed.BEDFileLoader;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.vcf.VCFFileLoader;
import com.hartwig.hmftools.common.variant.vcf.VCFSomaticFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class LoadSomaticVariants {

    private static final Logger LOGGER = LogManager.getLogger(LoadSomaticVariants.class);

    private static final String VCF_FILE = "vcf_file";
    private static final String REF_GENOME = "ref_genome";
    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    public static void main(@NotNull final String[] args) throws ParseException, IOException, HartwigException, SQLException {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String vcfFileLocation = cmd.getOptionValue(VCF_FILE);
        final String bedFileLocation = cmd.getOptionValue(HIGH_CONFIDENCE_BED);
        final String fastaFileLocation = cmd.getOptionValue(REF_GENOME);
        final DatabaseAccess dbAccess = databaseAccess(cmd);

        LOGGER.info("Reading somatic VCF File");
        final VCFSomaticFile vcfFile = VCFFileLoader.loadSomaticVCF(vcfFileLocation);
        final String sample = vcfFile.sample();

        LOGGER.info("Reading high confidence bed file");
        final Multimap<String, GenomeRegion> highConfidenceRegions = BEDFileLoader.fromBedFile(bedFileLocation);

        LOGGER.info("Loading indexed fasta reference file");
        IndexedFastaSequenceFile indexedFastaSequenceFile = new IndexedFastaSequenceFile(new File(fastaFileLocation));

        LOGGER.info("Querying purple database");
        final FittedPurity fittedPurity = dbAccess.readFittedPurity(sample);

        if (fittedPurity == null) {
            LOGGER.warn("Unable to retrieve purple data. Enrichment may be incomplete.");
        }
        final double purity = fittedPurity == null ? 1 : fittedPurity.purity();
        final double normFactor = fittedPurity == null ? 1 : fittedPurity.normFactor();

        final Multimap<String, PurpleCopyNumber> copyNumbers =
                Multimaps.index(dbAccess.readCopynumbers(sample), PurpleCopyNumber::chromosome);

        LOGGER.info("Enriching variants");
        final EnrichedSomaticVariantFactory enrichedSomaticVariantFactory =
                new EnrichedSomaticVariantFactory(purity, normFactor, highConfidenceRegions, copyNumbers, indexedFastaSequenceFile);
        final List<EnrichedSomaticVariant> variants = enrichedSomaticVariantFactory.enrich(vcfFile.variants());

        LOGGER.info("Persisting variants to database");
        dbAccess.writeSomaticVariants(sample, variants);

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createBasicOptions() {
        final Options options = new Options();
        options.addOption(REF_GENOME, true, "Path to the ref genome fasta file.");
        options.addOption(VCF_FILE, true, "Path to the vcf file.");
        options.addOption(HIGH_CONFIDENCE_BED, true, "Path to the high confidence bed file.");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private static DatabaseAccess databaseAccess(CommandLine cmd) throws SQLException {
        final String userName = cmd.getOptionValue(DB_USER);
        final String password = cmd.getOptionValue(DB_PASS);
        final String databaseUrl = cmd.getOptionValue(DB_URL);  //e.g. mysql://localhost:port/database";
        final String jdbcUrl = "jdbc:" + databaseUrl;
        return new DatabaseAccess(userName, password, jdbcUrl);
    }
}
