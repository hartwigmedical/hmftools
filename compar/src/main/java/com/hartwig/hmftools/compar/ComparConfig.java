package com.hartwig.hmftools.compar;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL;
import static com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_DESC;
import static com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig.addGenePanelOption;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.IGNORE_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadDelimitedIdFile;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.compar.common.CategoryType.ALL_CATEGORIES;
import static com.hartwig.hmftools.compar.common.CategoryType.DRIVER;
import static com.hartwig.hmftools.compar.common.CategoryType.GENE_COPY_NUMBER;
import static com.hartwig.hmftools.compar.common.CategoryType.LINX_CATEGORIES;
import static com.hartwig.hmftools.compar.common.CategoryType.PANEL_CATEGORIES;
import static com.hartwig.hmftools.compar.common.CategoryType.PURPLE_CATEGORIES;
import static com.hartwig.hmftools.compar.common.CategoryType.purpleCategories;
import static com.hartwig.hmftools.compar.common.CategoryType.linxCategories;
import static com.hartwig.hmftools.compar.common.FileSources.fromConfig;
import static com.hartwig.hmftools.compar.common.FileSources.registerConfig;
import static com.hartwig.hmftools.compar.common.MatchLevel.REPORTABLE;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.DB_DEFAULT_ARGS;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneFile;
import com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ComparConfig
{
    public final List<String> SampleIds;

    public final Map<CategoryType,MatchLevel> Categories;

    public final List<String> SourceNames; // list of sources to compare, eg prod vs pilot, or pipeline_1 vs pipeline_2

    public final Map<String,DatabaseAccess> DbConnections; // database access details keyed by source
    public final Map<String,FileSources> FileSources; // directories per type and keyed by source
    public final boolean RequiresLiftover;

    public final Set<String> DriverGenes;
    public final Set<String> IgnoreGenes;
    public final Set<String> AlternateTranscriptDriverGenes;
    public final boolean RestrictToDrivers;

    public final DiffThresholds Thresholds;

    public final String OutputDir;
    public final String OutputId;
    public final String ExpectedMismatchFile;

    public final boolean WriteDetailed;
    public final boolean IncludeMatches;
    public final int Threads;

    public final GenomeLiftoverCache LiftoverCache;

    private final Map<String,SampleIdMapping> mSampleIdMappings; // if required, mapping from original to new sampleId
    private boolean mIsValid;

    // config strings
    public static final String CATEGORIES = "categories";
    public static final String MATCH_LEVEL = "match_level";

    public static final String DB_SOURCE = "db_source";
    public static final String THRESHOLDS = "thresholds";

    public static final String WRITE_DETAILED_FILES = "write_detailed";
    public static final String INCLUDE_MATCHES = "include_matches";
    public static final String RESTRICT_TO_DRIVERS = "restrict_to_drivers";
    public static final String EXPECTED_MISMATCH_FILE = "expected_mismatches";
    public static final String IGNORE_GENES = "ignore_genes";

    public static final Logger CMP_LOGGER = LogManager.getLogger(ComparConfig.class);

    public static final String REF_SOURCE = "ref";
    public static final String NEW_SOURCE = "new";
    public static final String REQUIRES_LIFTOVER = "liftover";

    public ComparConfig(final ConfigBuilder configBuilder)
    {
        mIsValid = true;

        SampleIds = Lists.newArrayList();
        mSampleIdMappings = Maps.newHashMap();

        Categories = Maps.newHashMap();

        MatchLevel matchLevel = MatchLevel.valueOf(configBuilder.getValue(MATCH_LEVEL));

        String categoriesStr = configBuilder.getValue(CATEGORIES);

        CMP_LOGGER.info("default match level {}, categories: {}", matchLevel, categoriesStr);

        if(categoriesStr.equals(ALL_CATEGORIES))
        {
            Arrays.stream(CategoryType.values()).forEach(x -> Categories.put(x, matchLevel));
        }
        else if(categoriesStr.contains(PANEL_CATEGORIES))
        {
            if(categoriesStr.contains(PANEL_CATEGORIES))
                CategoryType.panelCategories().forEach(x -> Categories.put(x, matchLevel));
        }
        else
        {
            String itemDelim = categoriesStr.contains(ITEM_DELIM) ? ITEM_DELIM : CSV_DELIM;
            final String[] categoryStrings = categoriesStr.split(itemDelim);

            for(String catStr : categoryStrings)
            {
                if(catStr.equals(PURPLE_CATEGORIES))
                    purpleCategories().forEach(x -> Categories.put(x, matchLevel));
                else if(catStr.equals(LINX_CATEGORIES))
                    linxCategories().forEach(x -> Categories.put(x, matchLevel));
                else
                    Categories.put(CategoryType.valueOf(catStr), matchLevel);
            }
        }

        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);
        ExpectedMismatchFile = configBuilder.getValue(EXPECTED_MISMATCH_FILE);
        WriteDetailed = configBuilder.hasFlag(WRITE_DETAILED_FILES);
        IncludeMatches = configBuilder.hasFlag(INCLUDE_MATCHES);
        Threads = parseThreads(configBuilder);

        SourceNames = Lists.newArrayList(REF_SOURCE, NEW_SOURCE);
        loadSampleIds(configBuilder);

        RequiresLiftover = configBuilder.hasFlag(REQUIRES_LIFTOVER);

        DbConnections = Maps.newHashMap();
        FileSources = Maps.newHashMap();

        if(configBuilder.hasValue(formConfigSourceStr(DB_SOURCE, REF_SOURCE)) && configBuilder.hasValue(formConfigSourceStr(DB_SOURCE, NEW_SOURCE)))
        {
            loadDatabaseSources(configBuilder);
        }
        else
        {
            if(!loadFileSources(configBuilder))
            {
                mIsValid = false;
                CMP_LOGGER.error("missing DB or file source ref and new config");
            }
        }

        Thresholds = new DiffThresholds();
        Thresholds.loadConfig(configBuilder.getValue(THRESHOLDS, ""));

        DriverGenes = Sets.newHashSet();
        AlternateTranscriptDriverGenes = Sets.newHashSet();

        if(configBuilder.hasValue(DRIVER_GENE_PANEL))
        {
            try
            {
                List<DriverGene> driverGenes = DriverGeneFile.read(configBuilder.getValue(DRIVER_GENE_PANEL));

                for(DriverGene driverGene : driverGenes)
                {
                    DriverGenes.add(driverGene.gene());

                    if(!driverGene.additionalReportedTranscripts().isEmpty())
                        AlternateTranscriptDriverGenes.add(driverGene.gene());
                }
            }
            catch(IOException e)
            {
                CMP_LOGGER.error("failed to load driver gene panel file: {}", e.toString());
            }
        }

        IgnoreGenes = Sets.newHashSet();

        if(configBuilder.hasValue(IGNORE_GENES))
        {
            List<String> ignoredGenes = loadDelimitedIdFile(configBuilder.getValue(IGNORE_GENES), FLD_GENE_NAME, TSV_DELIM);
            ignoredGenes.forEach(x -> IgnoreGenes.add(x));
        }

        RestrictToDrivers = !DriverGenes.isEmpty() && configBuilder.hasFlag(RESTRICT_TO_DRIVERS);

        if(RestrictToDrivers)
        {
            CMP_LOGGER.info("restricting comparison to {} driver genes", DriverGenes.size());
        }

        LiftoverCache = new GenomeLiftoverCache(RequiresLiftover);
    }

    public boolean runCopyNumberGeneComparer() { return Categories.containsKey(DRIVER) && !Categories.containsKey(GENE_COPY_NUMBER); }

    public String sourceSampleId(final String source, final String sampleId)
    {
        if(mSampleIdMappings.isEmpty())
            return sampleId;

        SampleIdMapping mapping = mSampleIdMappings.get(sampleId);

        if(mapping == null || !mapping.SourceMapping.containsKey(source))
        {
            return sampleId;
        }

        return mapping.SourceMapping.get(source);
    }

    public String sourceGermlineSampleId(final String source, final String sampleId)
    {
        if(!mSampleIdMappings.isEmpty())
        {
            SampleIdMapping mapping = mSampleIdMappings.get(sampleId);

            if(mapping != null && mapping.GermlineSourceMapping.containsKey(source))
            {
                return mapping.GermlineSourceMapping.get(source);
            }
            else if(mapping != null && mapping.GermlineSampleId != null)
            {
                return mapping.GermlineSampleId;
            }
        }
        return sourceSampleId(source, sampleId) + "-ref";
    }

    public boolean isValid() { return mIsValid; }
    public boolean singleSample() { return SampleIds.size() == 1; }
    public boolean multiSample() { return SampleIds.size() > 1; }

    private static class SampleIdMapping
    {
        public final String SampleId;
        public final String GermlineSampleId;
        public Map<String,String> SourceMapping;
        public Map<String,String> GermlineSourceMapping;

        public SampleIdMapping(final String sampleId, final String germlineSampleId)
        {
            SampleId = sampleId;
            GermlineSampleId = germlineSampleId;
            SourceMapping = Maps.newHashMap();
            GermlineSourceMapping = Maps.newHashMap();
        }
    }

    private static final String COL_SAMPLE_ID = "SampleId";
    private static final String COL_GERMLINE_SAMPLE_ID = "GermlineSampleId";
    private static final String COL_REF_SAMPLE_ID = "RefSampleId";
    private static final String COL_REF_GERMLINE_SAMPLE_ID = "RefGermlineSampleId";
    private static final String COL_NEW_SAMPLE_ID = "NewSampleId";
    private static final String COL_NEW_GERMLINE_SAMPLE_ID = "NewGermlineSampleId";

    private void loadSampleIds(final ConfigBuilder configBuilder)
    {
        if(configBuilder.hasValue(SAMPLE_ID_FILE) && (configBuilder.hasValue(SAMPLE) || configBuilder.hasValue(REFERENCE)))
        {
            CMP_LOGGER.error("when the argument '{}' is set, the arguments '{}' and '{}' should not be set",
                    SAMPLE_ID_FILE, SAMPLE, REFERENCE);
            mIsValid = false;
            return;
        }

        if(configBuilder.hasValue(SAMPLE))
        {
            registerSampleIds(configBuilder.getValue(SAMPLE), configBuilder.getValue(REFERENCE, null));
            return;
        }

        if(!configBuilder.hasValue(SAMPLE_ID_FILE))
        {
            CMP_LOGGER.error("missing sample_id_file or sample config");
            mIsValid = false;
            return;
        }

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(configBuilder.getValue(SAMPLE_ID_FILE)));
            String header = lines.get(0);
            lines.remove(0);

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, CSV_DELIM);

            int sampleIndex = fieldsIndexMap.get(COL_SAMPLE_ID);
            Integer germlineSampleIndex = fieldsIndexMap.get(COL_GERMLINE_SAMPLE_ID);
            Integer refSampleIndex = fieldsIndexMap.get(COL_REF_SAMPLE_ID);
            Integer refGermlineSampleIndex = fieldsIndexMap.get(COL_REF_GERMLINE_SAMPLE_ID);
            Integer newSampleIndex = fieldsIndexMap.get(COL_NEW_SAMPLE_ID);
            Integer newGermlineSampleIndex = fieldsIndexMap.get(COL_NEW_GERMLINE_SAMPLE_ID);

            for(String line : lines)
            {
                if(line.isEmpty() || line.startsWith(IGNORE_SAMPLE_ID))
                    continue;

                String[] values = line.split(CSV_DELIM, -1);

                String sampleId = values[sampleIndex];
                String germlineSampleId = germlineSampleIndex != null ? values[germlineSampleIndex] : null;
                String refSampleId = refSampleIndex != null ? values[refSampleIndex] : null;
                String refGermlineSampleId = refGermlineSampleIndex != null ? values[refGermlineSampleIndex] : null;
                String newSampleId = newSampleIndex != null ? values[newSampleIndex] : null;
                String newGermlineSampleId = newGermlineSampleIndex != null ? values[newGermlineSampleIndex] : null;

                registerSampleIds(sampleId, germlineSampleId, refSampleId, refGermlineSampleId, newSampleId, newGermlineSampleId);
            }

            CMP_LOGGER.info("loaded {} samples from file", SampleIds.size());
        }
        catch(IOException e)
        {
            CMP_LOGGER.error("failed to load sample IDs: {}", e.toString());
        }
    }

    private void registerSampleIds(final String sampleId, final String germlineSampleId)
    {
        registerSampleIds(sampleId, germlineSampleId, null, null, null, null);
    }
    
    private void registerSampleIds(
            final String sampleId, final String germlineSampleId, final String refSampleId,
            final String refGermlineSampleId, final String newSampleId, final String newGermlineSampleId)
    {
        SampleIds.add(sampleId);

        SampleIdMapping mapping = new SampleIdMapping(sampleId, germlineSampleId);
        mSampleIdMappings.put(sampleId, mapping);

        if(refSampleId != null && SourceNames.size() >= 1)
            mapping.SourceMapping.put(SourceNames.get(0), refSampleId);
        if(newSampleId != null && SourceNames.size() >= 2)
            mapping.SourceMapping.put(SourceNames.get(1), newSampleId);

        if(refGermlineSampleId != null && SourceNames.size() >= 1)
            mapping.GermlineSourceMapping.put(SourceNames.get(0), refGermlineSampleId);
        if(newGermlineSampleId != null && SourceNames.size() >= 2)
            mapping.GermlineSourceMapping.put(SourceNames.get(1), newGermlineSampleId);
    }

    private static String formConfigSourceStr(final String sourceType, final String sourceName)
    {
        return format("%s_%s", sourceType, sourceName);
    }

    private void loadDatabaseSources(final ConfigBuilder configBuilder)
    {
        if(!configBuilder.hasValue(formConfigSourceStr(DB_SOURCE, REF_SOURCE)) || !configBuilder.hasValue(formConfigSourceStr(DB_SOURCE, NEW_SOURCE)))
            return;

        // form DB1;db_url;db_user;db_pass DB2;db_url;db_user;db_pass etc

        for(String sourceName : SourceNames)
        {
            String dbConfigValue =  configBuilder.getValue(formConfigSourceStr(DB_SOURCE, sourceName));
            String[] dbItems = dbConfigValue.split(CSV_DELIM, -1);

            if(dbItems.length != 3)
            {
                CMP_LOGGER.error("invalid DB source config({})", dbConfigValue);
                mIsValid = false;
                return;
            }

            String dbUrl = "jdbc:" + dbItems[0] + DB_DEFAULT_ARGS;
            String dbUsername = dbItems[1];
            String dbPass = dbItems[2];

            try
            {
                final DatabaseAccess dbAccess = new DatabaseAccess(dbUsername, dbPass, dbUrl);
                DbConnections.put(sourceName, dbAccess);
            }
            catch(SQLException e)
            {
                mIsValid = false;
            }
        }
    }

    private boolean loadFileSources(final ConfigBuilder configBuilder)
    {
        for(String sourceName : SourceNames)
        {
            FileSources fileSources = fromConfig(sourceName, configBuilder);

            if(fileSources == null)
            {
                FileSources.clear();
                return false;
            }

            FileSources.put(fileSources.Source, fileSources);
        }

        return true;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(
                CATEGORIES, false,
                "Categories to check separated by ';' from: PURITY, DRIVER, SOMATIC_VARIANT, GERMLINE_VARIANT, FUSION, DISRUPTION",
                ALL_CATEGORIES);

        configBuilder.addConfigItem(
                MATCH_LEVEL, false, "Match level from REPORTABLE (default) or DETAILED", REPORTABLE.toString());

        configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC);
        configBuilder.addConfigItem(REFERENCE, false, REFERENCE_DESC);
        addSampleIdFile(configBuilder, false);
        addGenePanelOption(configBuilder, false);
        configBuilder.addPath(IGNORE_GENES, false, "Genes to ignore in all comparisons, file with 'GeneName'");
        configBuilder.addConfigItem(THRESHOLDS, "In form: Field,AbsoluteDiff,PercentDiff, separated by ';'");

        configBuilder.addConfigItem(formConfigSourceStr(DB_SOURCE, REF_SOURCE), false, "Database configurations for reference data");
        configBuilder.addConfigItem(formConfigSourceStr(DB_SOURCE, NEW_SOURCE), false, "Database configurations for new data");

        registerConfig(configBuilder);

        configBuilder.addFlag(WRITE_DETAILED_FILES, "Write per-type details files");
        configBuilder.addConfigItem(EXPECTED_MISMATCH_FILE, "Existing expected mismatch file");
        configBuilder.addFlag(INCLUDE_MATCHES, "Also write matches to output file(s)");
        configBuilder.addFlag(RESTRICT_TO_DRIVERS, "Restrict any comparison involving genes to driver gene panel");
        configBuilder.addFlag(REQUIRES_LIFTOVER, "Lift over ref positions from v37 to v 38");

        addDatabaseCmdLineArgs(configBuilder, false);
        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
    }

    public ComparConfig()
    {
        mIsValid = true;

        SampleIds = Lists.newArrayList();
        Categories = Maps.newHashMap();
        OutputDir = null;
        OutputId = "";
        WriteDetailed = false;
        IncludeMatches = false;
        Threads = 0;

        DbConnections = Maps.newHashMap();
        FileSources = Maps.newHashMap();
        SourceNames = Lists.newArrayList();

        Thresholds = new DiffThresholds();
        DriverGenes = Sets.newHashSet();
        IgnoreGenes = Collections.emptySet();
        AlternateTranscriptDriverGenes = Sets.newHashSet();
        RestrictToDrivers = false;
        mSampleIdMappings = Maps.newHashMap();
        LiftoverCache = new GenomeLiftoverCache();
        RequiresLiftover = false;
        ExpectedMismatchFile = null;
    }
}
