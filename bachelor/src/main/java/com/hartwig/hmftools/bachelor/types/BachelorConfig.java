package com.hartwig.hmftools.bachelor.types;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.bachelor.GermlineVcfParser;
import com.hartwig.hmftools.bachelor.datamodel.Program;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class BachelorConfig {
    public final String SampleId;
    public final String OutputDir;
    public final String GermlineVcf;
    public final String BamFile;
    public final String RefGenomeFile;
    public final String PurpleDataDir;
    public final boolean SkipEnrichment;
    public final boolean IncludeVcfFiltered;

    public final Map<String, Program> ProgramConfigMap;

    private boolean mIsValid;

    // config options
    public static final String CONFIG_XML = "xml_config";

    // XML config regeneration:
    // xjc -d ./bachelor/target/generated-sources/xjc -p nl.hartwigmedicalfoundation.bachelor ./bachelor/src/main/resources/bachelor.xsd

    public static final String SAMPLE = "sample";

    public static final String DB_USER = "db_user";
    public static final String DB_PASS = "db_pass";
    public static final String DB_URL = "db_url";
    public static final String REF_GENOME = "ref_genome";

    private static final String GERMLINE_VCF = "germline_vcf";
    private static final String TUMOR_BAM_FILE = "tumor_bam_file";
    private static final String OUTPUT_DIR = "output_dir";
    private static final String PURPLE_DATA_DIRECTORY = "purple_data_dir";
    private static final String SKIP_ENRICHMENT = "skip_enrichment";
    private static final String INCLUDE_VCF_FILTERED = "include_vcf_filtered";

    public static final String LOG_DEBUG = "log_debug";

    public static final Logger BACH_LOGGER = LogManager.getLogger(BachelorConfig.class);

    public BachelorConfig(final CommandLine cmd)
    {
        mIsValid = true;

        ProgramConfigMap = Maps.newHashMap();

        if (cmd.hasOption(CONFIG_XML))
        {
            if (!loadXML(Paths.get(cmd.getOptionValue(CONFIG_XML)), ProgramConfigMap))
            {
                mIsValid = false;
            }
        }

        SampleId = cmd.getOptionValue(SAMPLE, "");

        if (SampleId.isEmpty()) {
            BACH_LOGGER.error("Missing Sample config");
            mIsValid = false;
        }

        GermlineVcf = cmd.getOptionValue(GERMLINE_VCF);

        SkipEnrichment = cmd.hasOption(SKIP_ENRICHMENT);
        IncludeVcfFiltered = cmd.hasOption(INCLUDE_VCF_FILTERED);
        BamFile = cmd.getOptionValue(TUMOR_BAM_FILE);
        RefGenomeFile = cmd.getOptionValue(REF_GENOME);

        String sampleOutputDir = cmd.getOptionValue(OUTPUT_DIR);

        if (!sampleOutputDir.endsWith(File.separator))
        {
            sampleOutputDir += File.separator;
        }

        OutputDir = sampleOutputDir;

        PurpleDataDir = cmd.getOptionValue(PURPLE_DATA_DIRECTORY, "");

        if (GermlineVcf == null)
        {
            BACH_LOGGER.error("Missing germline VCF file");
            mIsValid = false;
        }

        if (RefGenomeFile == null)
        {
            BACH_LOGGER.error("Missing ref genome input file");
            mIsValid = false;
        }

        if (!SkipEnrichment && (BamFile == null || PurpleDataDir.isEmpty()))
        {
            BACH_LOGGER.error("missing input files: BAM({}) purpleDataDir({})", BamFile, PurpleDataDir);

            mIsValid = false;
        }
    }

    public boolean isValid() {
        return mIsValid;
    }

    public static boolean loadXML(final Path path, Map<String, Program> configMap)
    {
        try
        {
            final ConfigSchema schema = ConfigSchema.make();

            final List<Program> programs = Files.walk(path)
                    .filter(p -> p.toString().endsWith(".xml"))
                    .map(schema::processXML)
                    .filter(Objects::nonNull)
                    .collect(Collectors.toList());

            for (final Program p : programs)
            {
                if (configMap.containsKey(p.getName()))
                {
                    BACH_LOGGER.error("Duplicate Programs detected: {}", p.getName());
                    return false;
                }
                else
                {
                    configMap.put(p.getName(), p);
                }
            }
        }
        catch (Exception e)
        {
            BACH_LOGGER.error("Error loading XML: {}", e.toString());
            return false;
        }

        return true;
    }

    @NotNull
    public static Options createOptions()
    {
        final Options options = new Options();

        // germline VCF parsing
        options.addOption(CONFIG_XML, true, "XML with genes, black and white lists");
        options.addOption(OUTPUT_DIR, true, "When in single-sample mode, all output written to this dir");
        options.addOption(TUMOR_BAM_FILE, true, "Location of a specific BAM file");
        options.addOption(GERMLINE_VCF, true, "Germline VCF file");
        options.addOption(SAMPLE, true, "Sample Id (not applicable for batch mode)");
        options.addOption(REF_GENOME, true, "Path to the ref genome fasta file");
        options.addOption(PURPLE_DATA_DIRECTORY, true, "Sub-directory with sample path for purple data");
        options.addOption(SKIP_ENRICHMENT, false, "Skip tumor BAM depth and Purple enrichment");
        options.addOption(INCLUDE_VCF_FILTERED, false, "Include variants filtered in the VCF");

        options.addOption(DB_USER, true, "Database user name");
        options.addOption(DB_PASS, true, "Database password");
        options.addOption(DB_URL, true, "Database url");

        // logging
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");

        GermlineVcfParser.addCmdLineOptions(options);

        return options;
    }
    
    @Nullable
    public static DatabaseAccess databaseAccess(@NotNull final CommandLine cmd)
    {
        if (!cmd.hasOption(DB_URL))
            return null;

        try
        {
            final String userName = cmd.getOptionValue(DB_USER);
            final String password = cmd.getOptionValue(DB_PASS);
            final String databaseUrl = cmd.getOptionValue(DB_URL);
            final String jdbcUrl = "jdbc:" + databaseUrl;
            return new DatabaseAccess(userName, password, jdbcUrl);
        }
        catch (SQLException e)
        {
            BACH_LOGGER.error("DB connection failed: {}", e.toString());
            return null;
        }
    }
}
