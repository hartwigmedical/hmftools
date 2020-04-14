package com.hartwig.hmftools.patientdb;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.region.BEDFileLoader;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.enrich.CompoundEnrichment;
import com.hartwig.hmftools.common.variant.enrich.HotspotEnrichment;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichmentLoadSomatics;
import com.hartwig.hmftools.common.variant.filter.SomaticFilter;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;

public class LoadSomaticVariants {

    private static final Logger LOGGER = LogManager.getLogger(LoadSomaticVariants.class);

    private static final String SAMPLE = "sample";

    private static final String SOMATIC_VCF = "somatic_vcf";
    private static final String HOTSPOT_TSV = "hotspot_tsv";
    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";

    private static final String PASS_FILTER = "pass_filter";
    private static final String SOMATIC_FILTER = "somatic_filter";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    public static void main(@NotNull String[] args) throws ParseException, IOException, SQLException {
        Options options = createBasicOptions();
        CommandLine cmd = createCommandLine(args, options);
        String vcfFileLocation = cmd.getOptionValue(SOMATIC_VCF);
        String highConfidenceBed = cmd.getOptionValue(HIGH_CONFIDENCE_BED);
        String sample = cmd.getOptionValue(SAMPLE);
        DatabaseAccess dbAccess = databaseAccess(cmd);
        CompoundFilter filter = new CompoundFilter(true);
        if (cmd.hasOption(PASS_FILTER)) {
            filter.add(new PassingVariantFilter());
        }
        if (cmd.hasOption(SOMATIC_FILTER)) {
            filter.add(new SomaticFilter());
        }

        CompoundEnrichment compoundEnrichment = new CompoundEnrichment();
        if (cmd.hasOption(HOTSPOT_TSV)) {
            String hotspotFile = cmd.getOptionValue(HOTSPOT_TSV);
            LOGGER.info("Reading hotspot file: {}", hotspotFile);
            compoundEnrichment.add(HotspotEnrichment.fromHotspotsFile(hotspotFile));
        }

        LOGGER.info("Reading high confidence bed file: {}", highConfidenceBed);
        Multimap<String, GenomeRegion> highConfidenceRegions = BEDFileLoader.fromBedFile(highConfidenceBed);

        LOGGER.info("Reading somatic VCF file: {}", vcfFileLocation);
        List<SomaticVariant> variants = SomaticVariantFactory.filteredInstanceWithEnrichment(filter,
                compoundEnrichment,
                VariantContextEnrichmentLoadSomatics.factory(highConfidenceRegions)).fromVCFFile(sample, vcfFileLocation);

        LOGGER.info("Persisting variants to database");
        dbAccess.writeSomaticVariants(sample, variants);

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createBasicOptions() {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample.");
        options.addOption(SOMATIC_VCF, true, "Path to the somatic SNV/indel vcf file.");
        options.addOption(HOTSPOT_TSV, true, "Location of hotspot tsv file");
        options.addOption(HIGH_CONFIDENCE_BED, true, "Path to the high confidence bed file.");
        options.addOption(PASS_FILTER, false, "Only load unfiltered variants");
        options.addOption(SOMATIC_FILTER, false, "Only load variants flagged SOMATIC");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        return new DefaultParser().parse(options, args);
    }

    @NotNull
    private static DatabaseAccess databaseAccess(@NotNull CommandLine cmd) throws SQLException {
        String userName = cmd.getOptionValue(DB_USER);
        String password = cmd.getOptionValue(DB_PASS);
        String databaseUrl = cmd.getOptionValue(DB_URL);  //e.g. mysql://localhost:port/database";
        String jdbcUrl = "jdbc:" + databaseUrl;
        return new DatabaseAccess(userName, password, jdbcUrl);
    }
}
