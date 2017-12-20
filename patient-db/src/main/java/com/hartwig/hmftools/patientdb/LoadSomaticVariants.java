package com.hartwig.hmftools.patientdb;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.bed.BEDFileLoader;
import com.hartwig.hmftools.common.region.bed.BEDFileLookup;
import com.hartwig.hmftools.common.region.bed.BEDFileLookupFileImpl;
import com.hartwig.hmftools.common.region.bed.BEDFileLookupNoImpl;
import com.hartwig.hmftools.common.variant.ClonalityCutoffKernel;
import com.hartwig.hmftools.common.variant.ClonalityFactory;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariant;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.filter.SomaticFilter;
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
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;

public class LoadSomaticVariants {

    private static final Logger LOGGER = LogManager.getLogger(LoadSomaticVariants.class);

    private static final String SAMPLE = "sample";
    private static final String VCF_FILE = "vcf_file";
    private static final String REF_GENOME = "ref_genome";
    private static final String PASS_FILTER = "pass_filter";
    private static final String SOMATIC_FILTER = "somatic_filter";
    private static final String MAPPABILITY_BED = "mappability_bed";
    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    public static void main(@NotNull final String[] args) throws ParseException, IOException, HartwigException, SQLException {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String vcfFileLocation = cmd.getOptionValue(VCF_FILE);
        final String highConfidenceBed = cmd.getOptionValue(HIGH_CONFIDENCE_BED);
        final String fastaFileLocation = cmd.getOptionValue(REF_GENOME);
        final String sample = cmd.getOptionValue(SAMPLE);
        final DatabaseAccess dbAccess = databaseAccess(cmd);
        final CompoundFilter filter = new CompoundFilter(true);
        if (cmd.hasOption(PASS_FILTER)) {
            filter.add(new PassingVariantFilter());
        }
        if (cmd.hasOption(SOMATIC_FILTER)) {
            filter.add(new SomaticFilter());
        }

        try (final BEDFileLookup mappabilityLookup = cmd.hasOption(MAPPABILITY_BED) ? new BEDFileLookupFileImpl(cmd.getOptionValue(
                MAPPABILITY_BED)) : new BEDFileLookupNoImpl()) {

            LOGGER.info("Reading somatic VCF File");
            final List<SomaticVariant> variants = new SomaticVariantFactory(filter).fromVCFFile(sample, vcfFileLocation);

            LOGGER.info("Reading high confidence bed file");
            final Multimap<String, GenomeRegion> highConfidenceRegions = BEDFileLoader.fromBedFile(highConfidenceBed);

            LOGGER.info("Loading indexed fasta reference file");
            IndexedFastaSequenceFile indexedFastaSequenceFile = new IndexedFastaSequenceFile(new File(fastaFileLocation));

            LOGGER.info("Querying purple database");
            final PurityContext purityContext = dbAccess.readPurityContext(sample);

            if (purityContext == null) {
                LOGGER.warn("Unable to retrieve purple data. Enrichment may be incomplete.");
            }

            final PurityAdjuster purityAdjuster = purityContext == null
                    ? new PurityAdjuster(Gender.FEMALE, 1, 1)
                    : new PurityAdjuster(purityContext.gender(), purityContext.bestFit().purity(), purityContext.bestFit().normFactor());

            final Multimap<String, PurpleCopyNumber> copyNumbers =
                    Multimaps.index(dbAccess.readCopynumbers(sample), PurpleCopyNumber::chromosome);

            final Multimap<String, FittedRegion> copyNumberRegions =
                    Multimaps.index(dbAccess.readCopyNumberRegions(sample), FittedRegion::chromosome);

            LOGGER.info("Incorporating purple purity");
            final PurityAdjustedSomaticVariantFactory purityAdjustmentFactory =
                    new PurityAdjustedSomaticVariantFactory(purityAdjuster, copyNumbers, copyNumberRegions);
            final List<PurityAdjustedSomaticVariant> purityAdjustedVariants = purityAdjustmentFactory.create(variants);

            final double clonalPloidy = ClonalityCutoffKernel.clonalCutoff(purityAdjustedVariants);

            LOGGER.info("Enriching variants");
            final EnrichedSomaticVariantFactory enrichedSomaticVariantFactory = new EnrichedSomaticVariantFactory(highConfidenceRegions,
                    indexedFastaSequenceFile,
                    mappabilityLookup,
                    new ClonalityFactory(purityAdjuster, clonalPloidy));
            final List<EnrichedSomaticVariant> enrichedVariants = enrichedSomaticVariantFactory.enrich(purityAdjustedVariants);

            LOGGER.info("Persisting variants to database");
            dbAccess.writeSomaticVariants(sample, enrichedVariants);
        }

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createBasicOptions() {
        final Options options = new Options();
        options.addOption(REF_GENOME, true, "Path to the ref genome fasta file.");
        options.addOption(VCF_FILE, true, "Path to the vcf file.");
        options.addOption(HIGH_CONFIDENCE_BED, true, "Path to the high confidence bed file.");
        options.addOption(MAPPABILITY_BED, true, "Path to the mappability score bed file.");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        options.addOption(SAMPLE, true, "Tumor sample.");
        options.addOption(PASS_FILTER, false, "Only load unfiltered variants");
        options.addOption(SOMATIC_FILTER, false, "Only load variants flagged SOMATIC");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static DatabaseAccess databaseAccess(@NotNull final CommandLine cmd) throws SQLException {
        final String userName = cmd.getOptionValue(DB_USER);
        final String password = cmd.getOptionValue(DB_PASS);
        final String databaseUrl = cmd.getOptionValue(DB_URL);  //e.g. mysql://localhost:port/database";
        final String jdbcUrl = "jdbc:" + databaseUrl;
        return new DatabaseAccess(userName, password, jdbcUrl);
    }
}
