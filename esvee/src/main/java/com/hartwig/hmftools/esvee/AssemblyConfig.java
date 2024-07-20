package com.hartwig.hmftools.esvee;

import static com.hartwig.hmftools.common.bam.BamToolName.BAMTOOL_PATH;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_DEBUG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_DEBUG_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_BAM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_BAMS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_BAM_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAMS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.parseLogReadIds;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.CONFIG_FILE_DELIM;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.esvee.AssemblyConstants.DEFAULT_ASSEMBLY_REF_BASE_WRITE_MAX;
import static com.hartwig.hmftools.esvee.alignment.BwaAligner.loadAlignerLibrary;
import static com.hartwig.hmftools.esvee.assembly.output.WriteType.ALIGNMENT_DATA;
import static com.hartwig.hmftools.esvee.assembly.output.WriteType.BREAKEND;
import static com.hartwig.hmftools.esvee.assembly.output.WriteType.fromConfig;
import static com.hartwig.hmftools.esvee.common.FileCommon.REF_GENOME_IMAGE_EXTENSION;
import static com.hartwig.hmftools.esvee.assembly.output.WriteType.ASSEMBLY_READ;
import static com.hartwig.hmftools.esvee.common.FileCommon.formEsveeInputFilename;
import static com.hartwig.hmftools.esvee.common.FileCommon.formPrepInputFilename;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.PREP_JUNCTIONS_FILE_ID;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.BamToolName;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.alignment.AlignmentCache;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.output.WriteType;
import com.hartwig.hmftools.esvee.common.ReadIdTrimmer;
import com.hartwig.hmftools.esvee.utils.TruthsetAnnotation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;

public class AssemblyConfig
{
    public final List<String> TumorIds;
    public final List<String> ReferenceIds;

    public final List<String> TumorBams;
    public final List<String> ReferenceBams;

    public final List<String> JunctionFiles;

    public final RefGenomeVersion RefGenVersion;
    public final RefGenomeCoordinates RefGenomeCoords;
    public final String RefGenomeFile;
    public final RefGenomeInterface RefGenome;
    public final String RefGenomeImageFile;
    public final String DecoyGenome;

    public final ValidationStringency BamStringency;

    public final List<WriteType> WriteTypes;

    public final boolean ProcessDiscordant;
    public final boolean RunAlignment;

    public final String OutputDir;
    public final String OutputId;
    public final String BamToolPath;

    public final SpecificRegions SpecificChrRegions;
    public final List<Junction> SpecificJunctions;

    public final boolean PerfDebug;
    public final double PerfLogTime;
    private final List<String> mLogReadIds;
    private final boolean mCheckLogReadIds;

    public final int AssemblyRefBaseWriteMax;
    public final int PhaseProcessingLimit;

    public final int Threads;

    public final String TruthsetFile;
    public final String AlignmentFile;

    public final boolean ApplyRemotePhasingReadCheckThreshold;

    private static final String REF_GENOME_IMAGE = "ref_genome_image";
    private static final String DECOY_GENOME = "decoy_genome";
    public static final String BWA_LIB_PATH = "bwa_lib";
    public static final String JUNCTION_FILES = "junction_files";

    private static final String WRITE_TYPES = "write_types";
    private static final String PERF_LOG_TIME = "perf_log_time";

    private static final String PROCESS_DISCORDANT = "discordant_pairs";
    private static final String RUN_ALIGNMENT = "run_alignment";

    private static final String PHASE_PROCESSING_LIMIT = "phase_process_limit";
    private static final String LOG_PHASE_GROUP_LINKS = "phase_group_links";
    private static final String SPECIFIC_JUNCTIONS = "specific_junctions";
    private static final String ASSEMBLY_REF_BASE_WRITE_MAX = "asm_ref_base_write_max";

    private static final String REMOTE_PHASING_READ_CHECK_THRESHOLD = "remote_phase_read_check_threshold";

    public static final Logger SV_LOGGER = LogManager.getLogger(AssemblyConfig.class);

    public static ReadIdTrimmer READ_ID_TRIMMER;

    public AssemblyConfig(final ConfigBuilder configBuilder)
    {
        TumorIds = Arrays.stream(configBuilder.getValue(TUMOR).split(CONFIG_FILE_DELIM)).collect(Collectors.toList());
        TumorBams = Arrays.stream(configBuilder.getValue(TUMOR_BAM).split(CONFIG_FILE_DELIM)).collect(Collectors.toList());

        if(configBuilder.hasValue(REFERENCE) && configBuilder.hasValue(REFERENCE_BAM))
        {
            ReferenceIds = Arrays.stream(configBuilder.getValue(REFERENCE).split(CONFIG_FILE_DELIM)).collect(Collectors.toList());
            ReferenceBams = Arrays.stream(configBuilder.getValue(REFERENCE_BAM).split(CONFIG_FILE_DELIM)).collect(Collectors.toList());
        }
        else
        {
            ReferenceIds = Collections.emptyList();
            ReferenceBams = Collections.emptyList();
        }

        if(TumorIds.isEmpty() || TumorIds.size() != TumorBams.size() || ReferenceIds.size() != ReferenceBams.size())
        {
            SV_LOGGER.error("tumor and reference IDs must match BAM files");
            System.exit(1);
        }

        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);

        JunctionFiles = Lists.newArrayList();

        if(configBuilder.hasValue(JUNCTION_FILES))
        {
            Arrays.stream(configBuilder.getValue(JUNCTION_FILES).split(CONFIG_FILE_DELIM)).forEach(x -> JunctionFiles.add(x));
        }
        else
        {
            // since Prep now reads multiple BAMs, only the tumor-labelled junctions file needs to be loaded
            List<String> combinedSampleIds = Lists.newArrayList(TumorIds);
            combinedSampleIds.addAll(ReferenceIds);

            for(String sampleId : combinedSampleIds)
            {
                String junctionFile = formPrepInputFilename(OutputDir, sampleId, PREP_JUNCTIONS_FILE_ID, OutputId);

                if(Files.exists(Paths.get(junctionFile)))
                    JunctionFiles.add(junctionFile);
            }
        }

        BamToolPath = configBuilder.getValue(BAMTOOL_PATH);

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        RefGenome = loadRefGenome(RefGenomeFile);

        RefGenomeImageFile = configBuilder.hasValue(REF_GENOME_IMAGE) ?
                configBuilder.getValue(REF_GENOME_IMAGE) : RefGenomeFile + REF_GENOME_IMAGE_EXTENSION;

        DecoyGenome = configBuilder.getValue(DECOY_GENOME);

        WriteTypes = fromConfig(configBuilder.getValue(WRITE_TYPES));

        AlignmentFile = AlignmentCache.filename(configBuilder);
        RunAlignment = configBuilder.hasFlag(RUN_ALIGNMENT) || AlignmentFile != null
                || WriteTypes.contains(BREAKEND) ||  WriteTypes.contains(ALIGNMENT_DATA);

        String bwaLibPath = configBuilder.getValue(BWA_LIB_PATH);

        if(RunAlignment || DecoyGenome != null)
            loadAlignerLibrary(bwaLibPath);

        ProcessDiscordant = configBuilder.hasFlag(PROCESS_DISCORDANT);
        BamStringency = ValidationStringency.STRICT;

        RefGenomeCoords = RefGenVersion == V37 ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;

        SpecificChrRegions = SpecificRegions.from(configBuilder);

        SpecificJunctions = Lists.newArrayList();
        if(configBuilder.hasValue(SPECIFIC_JUNCTIONS))
        {
            String[] specificJunctionsStr = configBuilder.getValue(SPECIFIC_JUNCTIONS).split(ITEM_DELIM);

            for(String specificJuncStr : specificJunctionsStr)
            {
                Junction junction = Junction.fromConfigStr(specificJuncStr);

                if(junction != null)
                    SpecificJunctions.add(junction);
            }

            SV_LOGGER.debug("loaded {} specific junctions", SpecificJunctions.size());
            Collections.sort(SpecificJunctions);
        }

        if(WriteTypes.contains(ASSEMBLY_READ) && !SpecificChrRegions.hasFilters() && SpecificJunctions.isEmpty())
        {
            SV_LOGGER.warn("writing assembly reads to TSV without region filtering may result in large output files & impact performance");
        }

        mLogReadIds = parseLogReadIds(configBuilder);
        mCheckLogReadIds = !mLogReadIds.isEmpty();

        PerfDebug = configBuilder.hasFlag(PERF_DEBUG);
        PerfLogTime = configBuilder.getDecimal(PERF_LOG_TIME);

        PhaseProcessingLimit = configBuilder.getInteger(PHASE_PROCESSING_LIMIT);

        // limit the length of ref bases by config unless using filters
        AssemblyRefBaseWriteMax = SpecificChrRegions.hasFilters() || !SpecificJunctions.isEmpty()
                ? 0 : configBuilder.getInteger(ASSEMBLY_REF_BASE_WRITE_MAX);

        Threads = parseThreads(configBuilder);

        TruthsetFile = TruthsetAnnotation.filename(configBuilder);

        ApplyRemotePhasingReadCheckThreshold = configBuilder.hasFlag(REMOTE_PHASING_READ_CHECK_THRESHOLD);

        READ_ID_TRIMMER = new ReadIdTrimmer(true);
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

    public void logReadId(final SAMRecord record, final String caller)
    {
        if(mCheckLogReadIds)
            logReadId(record.getReadName(), caller);
    }

    private void logReadId(final String readId, final String caller)
    {
        // debugging only
        if(mLogReadIds.contains(readId))
            SV_LOGGER.debug("caller({}) specific readId({})", caller, readId);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(TUMOR, true, TUMOR_IDS_DESC);
        configBuilder.addConfigItem(TUMOR_BAM, true, TUMOR_BAMS_DESC);

        configBuilder.addConfigItem(REFERENCE, false, REFERENCE_IDS_DESC);
        configBuilder.addConfigItem(REFERENCE_BAM, false, REFERENCE_BAMS_DESC);

        configBuilder.addPaths(
                JUNCTION_FILES, false, "List of SvPrep junction files, separated by ',', default is to match by sample name");

        addRefGenomeConfig(configBuilder, true);
        configBuilder.addPath(REF_GENOME_IMAGE, false, REFERENCE_BAM_DESC);
        configBuilder.addPath(DECOY_GENOME, false, "Decoy genome image file");

        configBuilder.addFlag(PROCESS_DISCORDANT, "Proces discordant-only groups");
        configBuilder.addFlag(RUN_ALIGNMENT, "Run assembly alignment");
        configBuilder.addPath(BWA_LIB_PATH, false, "Path to BWA library");

        String writeTypes = Arrays.stream(WriteType.values()).map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM));
        configBuilder.addConfigItem(WRITE_TYPES, false, "Write types from list: " + writeTypes);

        configBuilder.addConfigItem(LOG_READ_IDS, false, LOG_READ_IDS_DESC);
        configBuilder.addConfigItem(
                SPECIFIC_JUNCTIONS, false,
                "Specific junctions: format: Chromosome:Position:Orientation:Type (I, D or none), separated by ';'");

        configBuilder.addFlag(PERF_DEBUG, PERF_DEBUG_DESC);
        configBuilder.addDecimal(PERF_LOG_TIME, "Log performance data for routine exceeding specified time (0 = disabled)", 0);
        configBuilder.addFlag(LOG_PHASE_GROUP_LINKS, "Log assembly links to build phase groups");

        configBuilder.addInteger(
                ASSEMBLY_REF_BASE_WRITE_MAX, "Cap assembly ref bases in TSV and VCF, use zero to write all",
                DEFAULT_ASSEMBLY_REF_BASE_WRITE_MAX);

        configBuilder.addInteger(
                PHASE_PROCESSING_LIMIT, "Exclude phase groups above this size from extension and phase sets", 0);

        configBuilder.addFlag(REMOTE_PHASING_READ_CHECK_THRESHOLD, "Apply remote phase building max read check threshold");

        TruthsetAnnotation.registerConfig(configBuilder);
        AlignmentCache.registerConfig(configBuilder);
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

        RefGenVersion = V38;
        RefGenomeCoords = null;
        RefGenomeFile = null;
        RefGenome = null;
        RefGenomeImageFile = null;
        DecoyGenome = null;

        ProcessDiscordant = true;
        RunAlignment = true;

        BamStringency = ValidationStringency.SILENT;

        WriteTypes = Collections.emptyList();

        OutputDir = null;
        OutputId = null;
        BamToolPath = null;

        SpecificChrRegions = new SpecificRegions();
        SpecificJunctions = Collections.emptyList();

        PerfDebug = false;
        PerfLogTime = 0;
        mLogReadIds = Collections.emptyList();
        mCheckLogReadIds = false;

        AssemblyRefBaseWriteMax = 0;
        PhaseProcessingLimit = 0;
        Threads = 0;
        TruthsetFile = null;
        AlignmentFile = null;

        ApplyRemotePhasingReadCheckThreshold = false;

        READ_ID_TRIMMER = new ReadIdTrimmer(false);
    }
}
