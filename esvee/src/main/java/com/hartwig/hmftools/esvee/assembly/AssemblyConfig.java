package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.bamops.BamToolName.BAMTOOL_PATH;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_DEBUG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_DEBUG_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_BAM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_BAMS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAMS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.parseLogReadIds;
import static com.hartwig.hmftools.common.utils.config.ConfigItem.enumValueSelectionAsStr;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.pathFromFile;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.DEFAULT_ASSEMBLY_MAP_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.DEFAULT_ASSEMBLY_REF_BASE_WRITE_MAX;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.DEFAULT_DISC_RATE_INCREMENT;
import static com.hartwig.hmftools.esvee.assembly.SeqTechUtils.setSeqTechSpecifics;
import static com.hartwig.hmftools.esvee.common.FileCommon.JUNCTION_FILE;
import static com.hartwig.hmftools.esvee.common.FileCommon.JUNCTION_FILE_DESC;
import static com.hartwig.hmftools.esvee.common.FileCommon.PREP_DIR;
import static com.hartwig.hmftools.esvee.common.FileCommon.PREP_DIR_DESC;
import static com.hartwig.hmftools.esvee.common.FileCommon.REF_GENOME_IMAGE_EXTENSION;
import static com.hartwig.hmftools.esvee.common.FileCommon.formEsveeInputFilename;
import static com.hartwig.hmftools.esvee.common.FileCommon.formPrepBamFilenames;
import static com.hartwig.hmftools.esvee.common.FileCommon.formPrepInputFilename;
import static com.hartwig.hmftools.esvee.common.FileCommon.parseSampleBamLists;
import static com.hartwig.hmftools.esvee.common.FileCommon.registerCommonConfig;
import static com.hartwig.hmftools.esvee.common.FileCommon.setSequencingType;
import static com.hartwig.hmftools.esvee.common.WriteType.ASSEMBLY_READ;
import static com.hartwig.hmftools.esvee.common.WriteType.WRITE_TYPES;
import static com.hartwig.hmftools.esvee.common.WriteType.registerWriteTypes;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.PREP_JUNCTION_FILE_ID;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.common.ReadIdTrimmer;
import com.hartwig.hmftools.esvee.common.WriteType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class AssemblyConfig
{
    public final List<String> TumorIds;
    public final List<String> ReferenceIds;

    public final List<String> TumorBams;
    public final List<String> ReferenceBams;

    public final List<String> JunctionFiles;
    public final String PrepDir;

    public final RefGenomeVersion RefGenVersion;
    public final RefGenomeCoordinates RefGenomeCoords;
    public final String RefGenomeFile;
    public final RefGenomeInterface RefGenome;
    public final String RefGenomeImageFile;
    public final String DecoyGenome;

    public final Set<WriteType> WriteTypes;

    public final String OutputDir;
    public final String OutputId;
    public final String BamToolPath;

    public final int Threads;

    // default value overrides
    public static int PhaseProcessingLimit;
    public final int AssemblyMapQualThreshold;
    public final boolean DiscordantOnlyDisabled;
    public final double DiscordantRateIncrement;

    // debug options
    public final SpecificRegions SpecificChrRegions;
    public final List<Junction> SpecificJunctions;

    private final List<String> mLogReadIds;

    public final int AssemblyRefBaseWriteMax;
    public final boolean AssemblyDetailedTsv;

    public static boolean WriteCandidateReads;
    public static boolean AssemblyBuildDebug = false;
    public static boolean AssemblyLogAllReads = false;
    public static boolean DevDebug = false;
    public static boolean PerfDebug;
    public static double PerfLogTime;

    public final boolean ApplyRemotePhasingReadCheckThreshold;

    private static final String DECOY_GENOME = "decoy_genome";

    private static final String PERF_LOG_TIME = "perf_log_time";

    private static final String PHASE_PROCESSING_LIMIT = "phase_process_limit";
    private static final String DISC_RATE_INCREMENT = "disc_rate_increment";
    private static final String SPECIFIC_JUNCTIONS = "specific_junctions";
    private static final String ASSEMBLY_MAP_QUAL_THRESHOLD = "asm_map_qual_threshold";
    private static final String ASSEMBLY_REF_BASE_WRITE_MAX = "asm_ref_base_write_max";
    private static final String ASSEMBLY_BUILD_DEBUG = "asm_build_debug";
    private static final String ASSEMBLY_LOG_ALL_READS = "asm_log_all_reads";
    private static final String DISC_ONLY_DISABLED = "disc_only_disabled";
    private static final String WRITE_CANDIDATE_READS = "write_candidate_reads";
    private static final String ASSEMBLY_TSV_DETAILED = "assembly_detailed_tsv";

    private static final String REMOTE_PHASING_READ_CHECK_THRESHOLD = "remote_phase_read_check_threshold";

    public static final Logger SV_LOGGER = LogManager.getLogger(AssemblyConfig.class);

    public static ReadIdTrimmer READ_ID_TRIMMER;

    public AssemblyConfig(final ConfigBuilder configBuilder)
    {
        this(configBuilder, false);
    }

    public AssemblyConfig(final ConfigBuilder configBuilder, boolean asSubRoutine)
    {
        String prepBamDir = null;

        if(!asSubRoutine)
        {
            if(configBuilder.hasValue(TUMOR_BAM))
            {
                List<String> tumorBams = parseSampleBamLists(configBuilder, TUMOR_BAM);
                prepBamDir = pathFromFile(tumorBams.get(0));
            }
            else if(configBuilder.hasValue(REFERENCE_BAM))
            {
                List<String> refBams = parseSampleBamLists(configBuilder, REFERENCE_BAM);
                prepBamDir = pathFromFile(refBams.get(0));
            }
        }

        if(!configBuilder.hasValue(OUTPUT_DIR) && prepBamDir != null)
        {
            OutputDir = prepBamDir;
        }
        else
        {
            OutputDir = parseOutputDir(configBuilder);
        }

        if(configBuilder.hasValue(PREP_DIR))
            PrepDir = checkAddDirSeparator(configBuilder.getValue(PREP_DIR));
        else if(prepBamDir != null)
            PrepDir = prepBamDir;
        else
            PrepDir = OutputDir;

        OutputId = configBuilder.getValue(OUTPUT_ID);

        TumorIds = parseSampleBamLists(configBuilder, TUMOR);

        TumorBams = Lists.newArrayList();
        ReferenceBams = Lists.newArrayList();
        ReferenceIds = Lists.newArrayList();

        List<String> prepTumorBams = formPrepBamFilenames(PrepDir, TumorIds);

        if(prepTumorBams.size() == TumorIds.size())
            TumorBams.addAll(prepTumorBams);
        else if(!asSubRoutine)
            TumorBams.addAll(parseSampleBamLists(configBuilder, TUMOR_BAM));

        if(configBuilder.hasValue(REFERENCE))
        {
            ReferenceIds.addAll(parseSampleBamLists(configBuilder, REFERENCE));

            List<String> prepRefBams = formPrepBamFilenames(PrepDir, ReferenceIds);

            if(prepRefBams.size() == ReferenceIds.size())
                ReferenceBams.addAll(prepRefBams);
            else if(!asSubRoutine)
                ReferenceBams.addAll(parseSampleBamLists(configBuilder, REFERENCE_BAM));
        }

        if(TumorIds.isEmpty() && ReferenceIds.isEmpty())
        {
            SV_LOGGER.error("no tumor or reference IDs provided");
            System.exit(1);
        }

        if(TumorIds.size() != TumorBams.size() || ReferenceIds.size() != ReferenceBams.size())
        {
            SV_LOGGER.error("tumor and reference IDs must match BAM files");
            System.exit(1);
        }

        if(asSubRoutine)
        {
            if(!TumorBams.isEmpty())
                SV_LOGGER.debug("processing tumor prep bam(s): {}", TumorBams);

            if(!ReferenceBams.isEmpty())
                SV_LOGGER.debug("processing reference prep bam(s): {}", ReferenceBams);
        }

        // in germline mode, transfer to reference values to the tumor ones
        if(TumorIds.isEmpty() && !ReferenceIds.isEmpty())
        {
            TumorIds.addAll(ReferenceIds);
            ReferenceIds.clear();

            TumorBams.addAll(ReferenceBams);
            ReferenceBams.clear();
        }

        JunctionFiles = Lists.newArrayList();

        if(configBuilder.hasValue(JUNCTION_FILE))
        {
            JunctionFiles.add(configBuilder.getValue(JUNCTION_FILE));
        }
        else
        {
            // since Prep now reads multiple BAMs, only the tumor-labelled junctions file needs to be loaded
            String junctionFile = formPrepInputFilename(PrepDir, TumorIds.get(0), PREP_JUNCTION_FILE_ID, OutputId);

            if(Files.exists(Paths.get(junctionFile)))
                JunctionFiles.add(junctionFile);
        }

        BamToolPath = configBuilder.getValue(BAMTOOL_PATH);

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        RefGenome = loadRefGenome(RefGenomeFile);
        RefGenomeImageFile = RefGenomeFile + REF_GENOME_IMAGE_EXTENSION;

        DecoyGenome = configBuilder.getValue(DECOY_GENOME);

        WriteTypes = WriteType.parseAssemblyTypes(configBuilder.getValue(WRITE_TYPES));

        setSequencingType(configBuilder);
        setSeqTechSpecifics();

        RefGenomeCoords = RefGenVersion == V37 ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;

        AssemblyMapQualThreshold = configBuilder.getInteger(ASSEMBLY_MAP_QUAL_THRESHOLD);

        SpecificChrRegions = SpecificRegions.from(configBuilder);

        SpecificJunctions = Lists.newArrayList();
        if(configBuilder.hasValue(SPECIFIC_JUNCTIONS))
        {
            String specificJunctionsStr = configBuilder.getValue(SPECIFIC_JUNCTIONS);
            String junctionDelim = specificJunctionsStr.contains("_") ? "_" : ITEM_DELIM;
            String[] specificJunctionsList = specificJunctionsStr.split(junctionDelim);

            for(String specificJuncStr : specificJunctionsList)
            {
                Junction junction = Junction.fromConfigStr(specificJuncStr);

                if(junction != null)
                    SpecificJunctions.add(junction);
            }

            SV_LOGGER.debug("loaded {} specific junctions", SpecificJunctions.size());
            Collections.sort(SpecificJunctions);
        }

        boolean hasFilters = SpecificChrRegions.hasFilters() || !SpecificJunctions.isEmpty();

        if(WriteTypes.contains(ASSEMBLY_READ) && !hasFilters)
        {
            SV_LOGGER.warn("writing assembly reads to TSV without region filtering may result in large output files & impact performance");
        }

        AssemblyDetailedTsv = configBuilder.hasFlag(ASSEMBLY_TSV_DETAILED) || hasFilters;

        mLogReadIds = parseLogReadIds(configBuilder);

        DiscordantOnlyDisabled = configBuilder.hasFlag(DISC_ONLY_DISABLED);

        DiscordantRateIncrement = configBuilder.getDecimal(DISC_RATE_INCREMENT);

        PerfLogTime = configBuilder.getDecimal(PERF_LOG_TIME);
        PerfDebug = configBuilder.hasFlag(PERF_DEBUG) || PerfLogTime > 0;
        AssemblyBuildDebug = configBuilder.hasFlag(ASSEMBLY_BUILD_DEBUG);
        AssemblyLogAllReads = configBuilder.hasFlag(ASSEMBLY_BUILD_DEBUG);
        WriteCandidateReads = configBuilder.hasFlag(WRITE_CANDIDATE_READS);

        PhaseProcessingLimit = configBuilder.getInteger(PHASE_PROCESSING_LIMIT);

        // limit the length of ref bases by config unless using filters
        AssemblyRefBaseWriteMax = hasFilters ? 0 : configBuilder.getInteger(ASSEMBLY_REF_BASE_WRITE_MAX);

        Threads = parseThreads(configBuilder);

        ApplyRemotePhasingReadCheckThreshold = configBuilder.hasFlag(REMOTE_PHASING_READ_CHECK_THRESHOLD);

        READ_ID_TRIMMER = new ReadIdTrimmer(!hasFilters && TumorBams.size() == 1);
    }

    public List<String> combinedSampleIds()
    {
        List<String> combinedSampleIds = Lists.newArrayList(TumorIds);
        combinedSampleIds.addAll(ReferenceIds);
        return combinedSampleIds;
    }

    public List<String> combinedBamFiles()
    {
        List<String> combinedSampleIds = Lists.newArrayList(TumorBams);
        combinedSampleIds.addAll(ReferenceBams);
        return combinedSampleIds;
    }

    public String sampleId() { return TumorIds.get(0); }

    public String outputFilename(final WriteType writeType)
    {
        return formEsveeInputFilename(OutputDir, sampleId(), writeType.fileId(), OutputId);
    }

    private void logReadId(final String readId, final String caller)
    {
        // debugging only
        if(mLogReadIds.contains(readId))
            SV_LOGGER.debug("caller({}) specific readId({})", caller, readId);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(TUMOR, false, TUMOR_IDS_DESC);
        configBuilder.addConfigItem(TUMOR_BAM, false, TUMOR_BAMS_DESC);

        configBuilder.addConfigItem(REFERENCE, false, REFERENCE_IDS_DESC);
        configBuilder.addConfigItem(REFERENCE_BAM, false, REFERENCE_BAMS_DESC);

        configBuilder.addPaths(JUNCTION_FILE, false, JUNCTION_FILE_DESC);
        configBuilder.addPaths(PREP_DIR, false, PREP_DIR_DESC);

        registerCommonConfig(configBuilder);
        registerWriteTypes(configBuilder);

        addRefGenomeConfig(configBuilder, true);
        configBuilder.addPath(DECOY_GENOME, false, "Decoy genome image file");

        if(!configBuilder.isRegistered(WRITE_TYPES))
            configBuilder.addConfigItem(WRITE_TYPES, false, enumValueSelectionAsStr(WriteType.values(), "Write types"));

        configBuilder.addConfigItem(LOG_READ_IDS, false, LOG_READ_IDS_DESC);
        configBuilder.addConfigItem(
                SPECIFIC_JUNCTIONS, false,
                "Specific junctions: format: Chromosome:Position:Orientation:Type (I, D or none), separated by ';'");

        configBuilder.addFlag(PERF_DEBUG, PERF_DEBUG_DESC);
        configBuilder.addDecimal(PERF_LOG_TIME, "Log performance data for routine exceeding specified time (0 = disabled)", 0);

        configBuilder.addInteger(
                ASSEMBLY_REF_BASE_WRITE_MAX, "Cap assembly ref bases in TSV and VCF, use zero to write all",
                DEFAULT_ASSEMBLY_REF_BASE_WRITE_MAX);

        configBuilder.addInteger(
                PHASE_PROCESSING_LIMIT, "Exclude phase groups above this size from extension and phase sets", 0);

        configBuilder.addFlag(DISC_ONLY_DISABLED, "Disable discordant-only junctions");

        configBuilder.addInteger(
                ASSEMBLY_MAP_QUAL_THRESHOLD, "Realign and test assemblies with average map-qual below this threshold",
                DEFAULT_ASSEMBLY_MAP_QUAL_THRESHOLD);

        configBuilder.addFlag(REMOTE_PHASING_READ_CHECK_THRESHOLD, "Apply remote phase building max read check threshold");
        configBuilder.addFlag(WRITE_CANDIDATE_READS, "Write assembly candidate reads regardless of whether used");

        configBuilder.addFlag(ASSEMBLY_BUILD_DEBUG, "Log assembly building working");
        configBuilder.addFlag(ASSEMBLY_LOG_ALL_READS, "Log assembly all reads even if not used in support");
        configBuilder.addFlag(ASSEMBLY_TSV_DETAILED, "Log extra assembly info to TSV");

        configBuilder.addDecimal(DISC_RATE_INCREMENT, "Discordant rate increment", DEFAULT_DISC_RATE_INCREMENT);

        BamToolName.addConfig(configBuilder);

        SpecificRegions.addSpecificChromosomesRegionsConfig(configBuilder);

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
    }

    @VisibleForTesting
    public AssemblyConfig()
    {
        TumorIds = Collections.emptyList();
        ReferenceIds = Collections.emptyList();

        TumorBams = Collections.emptyList();
        ReferenceBams = Collections.emptyList();

        JunctionFiles = Collections.emptyList();
        PrepDir = null;

        RefGenVersion = V38;
        RefGenomeCoords = null;
        RefGenomeFile = null;
        RefGenome = null;
        RefGenomeImageFile = null;
        DecoyGenome = null;

        WriteTypes = Collections.emptySet();

        OutputDir = null;
        OutputId = null;
        BamToolPath = null;

        SpecificChrRegions = new SpecificRegions();
        SpecificJunctions = Collections.emptyList();

        PerfDebug = false;
        PerfLogTime = 0;
        AssemblyDetailedTsv = false;
        mLogReadIds = Collections.emptyList();

        AssemblyMapQualThreshold = -1;
        AssemblyRefBaseWriteMax = 0;
        PhaseProcessingLimit = 0;
        DiscordantOnlyDisabled = false;
        Threads = 0;
        DiscordantRateIncrement = DEFAULT_DISC_RATE_INCREMENT;

        ApplyRemotePhasingReadCheckThreshold = false;
        AssemblyBuildDebug = false;
        WriteCandidateReads = false;

        READ_ID_TRIMMER = new ReadIdTrimmer(false);
    }
}
