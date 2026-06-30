package com.hartwig.hmftools.compar;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL;
import static com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig.addGenePanelOption;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigItem.enumValueSelectionAsStr;
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
import static com.hartwig.hmftools.compar.common.SourceType.NEW;
import static com.hartwig.hmftools.compar.common.SourceType.OLD;
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
import com.hartwig.hmftools.compar.common.SourceData;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.WriteType;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ComparConfig
{
    public final List<String> SampleIds;
    public final Map<String,String> SampleToReferenceIds; // mapping of tumor to reference ID

    public final Map<CategoryType,MatchLevel> Categories;

    public final List<SourceData> Sources;

    public final boolean RequiresLiftover;

    public final Set<String> DriverGenes;
    public final Set<String> IgnoreGenes;
    public final Set<String> AlternateTranscriptDriverGenes;
    public final boolean RestrictToDrivers;

    public final DiffThresholds Thresholds;

    public final String OutputDir;
    public final String OutputId;
    public final String KnownMismatchFile;

    public final List<WriteType> WriteTypes;
    public final boolean IncludeMatches;
    public final int Threads;

    public final GenomeLiftoverCache LiftoverCache;

    private boolean mIsValid;

    // config strings
    public static final String CATEGORIES = "categories";
    public static final String MATCH_LEVEL = "match_level";

    public static final String DB_SOURCE = "db_source";
    public static final String THRESHOLDS = "thresholds";

    public static final String WRITE_TYPES = "write_types";
    public static final String WRITE_DETAILED_FILES = "write_detailed";
    public static final String INCLUDE_MATCHES = "include_matches";
    public static final String RESTRICT_TO_DRIVERS = "restrict_to_drivers";
    public static final String KNOWN_MISMATCH_FILE = "known_mismatches";
    public static final String IGNORE_GENES = "ignore_genes";

    public static final String OLD_SOURCE_CFG = OLD.configStr();
    public static final String NEW_SOURCE_CFG = NEW.configStr();

    public static final Logger CMP_LOGGER = LogManager.getLogger(ComparConfig.class);

    public static final String REQUIRES_LIFTOVER = "liftover";

    public ComparConfig(final ConfigBuilder configBuilder)
    {
        mIsValid = true;

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
        KnownMismatchFile = configBuilder.getValue(KNOWN_MISMATCH_FILE);

        WriteTypes = Lists.newArrayList();

        if(configBuilder.hasFlag(WRITE_DETAILED_FILES))
        {
            WriteTypes.add(WriteType.TYPE_SPECIFIC);
        }
        else
        {
            WriteTypes.addAll(WriteType.fromConfigStr(configBuilder.getValue(WRITE_TYPES)));
        }

        IncludeMatches = configBuilder.hasFlag(INCLUDE_MATCHES);
        Threads = parseThreads(configBuilder);

        // build old and new data sources
        Sources = Lists.newArrayListWithCapacity(2);

        RequiresLiftover = configBuilder.hasFlag(REQUIRES_LIFTOVER);

        if(configBuilder.hasValue(formConfigSourceStr(DB_SOURCE, OLD.configStr()))
        && configBuilder.hasValue(formConfigSourceStr(DB_SOURCE, NEW.configStr())))
        {
            loadDatabaseSources(configBuilder);
        }
        else
        {
            if(!loadFileSources(configBuilder))
            {
                mIsValid = false;
                CMP_LOGGER.error("missing DB or file source old and new config");
            }
        }

        SampleIds = Lists.newArrayList();
        SampleToReferenceIds = Maps.newHashMap();

        loadSampleIds(configBuilder);

        Thresholds = new DiffThresholds();

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

            CMP_LOGGER.info("ignoring {} genes", IgnoreGenes.size());
        }

        RestrictToDrivers = !DriverGenes.isEmpty() && configBuilder.hasFlag(RESTRICT_TO_DRIVERS);

        if(RestrictToDrivers)
        {
            CMP_LOGGER.info("restricting comparison to {} driver genes", DriverGenes.size());
        }

        LiftoverCache = new GenomeLiftoverCache(RequiresLiftover);
    }

    public SourceData getSourceData(final SourceType sourceType)
    {
        return Sources.stream().filter(x -> x.Type == sourceType).findFirst().orElse(null);
    }

    public String sourceSampleId(final SourceType sourceType, final String sampleId)
    {
        SourceData sourceData = getSourceData(sourceType);
        return sourceData.SampleIdMapping.getOrDefault(sampleId, sampleId);
    }

    public String sourceReferenceId(final SourceType sourceType, final String sampleId)
    {
        SourceData sourceData = getSourceData(sourceType);
        String referenceId = SampleToReferenceIds.get(sampleId);
        return sourceData.ReferenceSampleIdMapping.getOrDefault(sampleId, referenceId);
    }

    public boolean isValid() { return mIsValid; }
    public boolean singleSample() { return SampleIds.size() == 1; }
    public boolean multiSample() { return SampleIds.size() > 1; }

    private enum SampleIdFileColumn
    {
        SampleId,
        ReferenceId,
        OldSampleId,
        OldReferenceId,
        NewSampleId,
        NewReferenceId;
    }

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
            String sampleId = configBuilder.getValue(SAMPLE);
            String referenceId = configBuilder.getValue(REFERENCE, sampleId + "-ref"); // Hartwig default naming
            registerSampleIds(sampleId, referenceId);
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

            int sampleIdIndex = fieldsIndexMap.get(SampleIdFileColumn.SampleId.toString());
            Integer referenceIdIndex = fieldsIndexMap.get(SampleIdFileColumn.ReferenceId.toString());
            Integer oldSampleIndex = fieldsIndexMap.get(SampleIdFileColumn.OldSampleId.toString());
            Integer oldReferenceIdIndex = fieldsIndexMap.get(SampleIdFileColumn.OldReferenceId.toString());
            Integer newSampleIdIndex = fieldsIndexMap.get(SampleIdFileColumn.NewSampleId.toString());
            Integer newReferenceIdndex = fieldsIndexMap.get(SampleIdFileColumn.NewReferenceId.toString());

            for(String line : lines)
            {
                if(line.isEmpty() || line.startsWith(IGNORE_SAMPLE_ID))
                    continue;

                String[] values = line.split(CSV_DELIM, -1);

                String sampleId = values[sampleIdIndex];
                String referenceId = referenceIdIndex != null ? values[referenceIdIndex] : null;
                String oldSampleId = oldSampleIndex != null ? values[oldSampleIndex] : null;
                String oldReferenceSampleId = oldReferenceIdIndex != null ? values[oldReferenceIdIndex] : null;
                String newSampleId = newSampleIdIndex != null ? values[newSampleIdIndex] : null;
                String newReferenceSampleId = newReferenceIdndex != null ? values[newReferenceIdndex] : null;

                registerSampleIds(sampleId, referenceId, oldSampleId, oldReferenceSampleId, newSampleId, newReferenceSampleId);
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
            final String sampleId, final String referenceId, final String oldSampleId,
            final String oldReferenceSampleId, final String newSampleId, final String newReferenceSampleId)
    {
        SampleIds.add(sampleId);
        SampleToReferenceIds.put(sampleId, referenceId);

        SourceData oldSourceData = getSourceData(OLD);
        SourceData newSourceData = getSourceData(NEW);

        if(oldSampleId != null)
            oldSourceData.SampleIdMapping.put(sampleId, oldSampleId);

        if(oldReferenceSampleId != null)
            oldSourceData.ReferenceSampleIdMapping.put(sampleId, oldReferenceSampleId);

        if(newSampleId != null)
            newSourceData.SampleIdMapping.put(sampleId, newSampleId);

        if(newReferenceSampleId != null)
            newSourceData.ReferenceSampleIdMapping.put(sampleId, newReferenceSampleId);
    }

    private static String formConfigSourceStr(final String sourceType, final String sourceName)
    {
        return format("%s_%s", sourceType, sourceName);
    }

    private void loadDatabaseSources(final ConfigBuilder configBuilder)
    {
        // form DB1;db_url;db_user;db_pass DB2;db_url;db_user;db_pass etc

        for(SourceType sourceType : SourceType.values())
        {
            String dbConfigValue =  configBuilder.getValue(formConfigSourceStr(DB_SOURCE, sourceType.configStr()));
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
                DatabaseAccess dbAccess = new DatabaseAccess(dbUsername, dbPass, dbUrl);
                Sources.add(new SourceData(sourceType, dbAccess, null));
            }
            catch(SQLException e)
            {
                mIsValid = false;
            }
        }
    }

    private boolean loadFileSources(final ConfigBuilder configBuilder)
    {
        for(SourceType sourceType : SourceType.values())
        {
            FileSources fileSources = fromConfig(sourceType, configBuilder);
            Sources.add(new SourceData(sourceType, null, fileSources));
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

        configBuilder.addConfigItem(
                formConfigSourceStr(DB_SOURCE, OLD_SOURCE_CFG), false, "Database configurations for reference data");

        configBuilder.addConfigItem(
                formConfigSourceStr(DB_SOURCE, NEW_SOURCE_CFG), false, "Database configurations for new data");

        registerConfig(configBuilder);

        configBuilder.addConfigItem(WRITE_TYPES, enumValueSelectionAsStr(WriteType.values(), "Write types"));
        configBuilder.addFlag(WRITE_DETAILED_FILES, "Write per-type details files");
        configBuilder.addConfigItem(KNOWN_MISMATCH_FILE, "File with sample curations or expected mismatches");
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
        SampleToReferenceIds = Maps.newHashMap();
        Sources = Lists.newArrayList();
        Sources.add(new SourceData(OLD, null, null));
        Sources.add(new SourceData(NEW, null, null));

        Categories = Maps.newHashMap();
        OutputDir = null;
        OutputId = "";
        IncludeMatches = false;
        Threads = 0;
        WriteTypes = WriteType.DEFAULT_WRITE_TYPES;

        Thresholds = new DiffThresholds();
        DriverGenes = Sets.newHashSet();
        IgnoreGenes = Collections.emptySet();
        AlternateTranscriptDriverGenes = Sets.newHashSet();
        RestrictToDrivers = false;
        LiftoverCache = new GenomeLiftoverCache();
        RequiresLiftover = false;
        KnownMismatchFile = null;
    }
}
