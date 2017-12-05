package com.hartwig.hmftools.svannotation;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.cosmicfusions.COSMICGeneFusionModel;
import com.hartwig.hmftools.common.cosmicfusions.COSMICGeneFusions;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory;
import com.hartwig.hmftools.hmfslicer.HmfGenePanelSupplier;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalysis;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalyzer;
import com.hartwig.hmftools.svannotation.dao.StructuralVariantAnnotationDAO;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class LoadStructuralVariants {

    private static final Logger LOGGER = LogManager.getLogger(LoadStructuralVariants.class);

    private static final String SAMPLE = "sample";
    private static final String VCF_FILE = "vcf_file";
    private static final String ENSEMBL_DB = "ensembl_db";
    private static final String FUSION_CSV = "fusion_csv";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    public static void main(@NotNull final String[] args) throws ParseException, IOException, HartwigException, SQLException {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String vcfFileLocation = cmd.getOptionValue(VCF_FILE);
        final DatabaseAccess dbAccess = databaseAccess(cmd);
        final StructuralVariantFactory factory = new StructuralVariantFactory();
        final String tumorSample = cmd.getOptionValue(SAMPLE);

        LOGGER.info("Reading VCF File");
        try (final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(vcfFileLocation,
                new VCFCodec(), false)) {
            reader.iterator().forEach(factory::addVariantContext);
        }

        LOGGER.info("Persisting variants to database");
        dbAccess.writeStructuralVariants(tumorSample, factory.results());

        LOGGER.info("Reading back...");
        final List<StructuralVariant> structuralVariantList = dbAccess.readStructuralVariants(tumorSample);

        final VariantAnnotator annotator = MySQLAnnotator.make("jdbc:" + cmd.getOptionValue(ENSEMBL_DB));
        final COSMICGeneFusionModel cosmicGeneFusions = COSMICGeneFusions.readFromCSV(cmd.getOptionValue(FUSION_CSV));

        LOGGER.info("Annotating...");
        final StructuralVariantAnalyzer analyzer =
                new StructuralVariantAnalyzer(annotator, HmfGenePanelSupplier.hmfGeneList(), cosmicGeneFusions);
        final StructuralVariantAnalysis analysis = analyzer.run(structuralVariantList);

        final StructuralVariantAnnotationDAO annotationDAO = new StructuralVariantAnnotationDAO(dbAccess.getContext());

        LOGGER.info("Persisting annotations to database...");
        annotationDAO.write(analysis);

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createBasicOptions() {
        final Options options = new Options();
        options.addOption(VCF_FILE, true, "Path to the vcf file.");
        options.addOption(SAMPLE, true, "Tumor sample.");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        options.addOption(FUSION_CSV, true, "Path towards a CSV containing white-listed gene fusions.");
        options.addOption(ENSEMBL_DB, true, "Annotate structural variants using this Ensembl DB URI");
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
