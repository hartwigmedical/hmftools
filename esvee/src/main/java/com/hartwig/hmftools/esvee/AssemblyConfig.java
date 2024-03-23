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
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.pathFromFile;
import static com.hartwig.hmftools.esvee.AssemblyConstants.DEFAULT_ASSEMBLY_REF_BASE_WRITE_MAX;
import static com.hartwig.hmftools.esvee.AssemblyConstants.REF_GENOME_IMAGE_EXTENSION;
import static com.hartwig.hmftools.esvee.AssemblyConstants.SV_PREP_JUNCTIONS_FILE_ID;
import static com.hartwig.hmftools.esvee.common.CommonUtils.formOutputFile;
import static com.hartwig.hmftools.esvee.output.WriteType.ASSEMBLY_BAM;
import static com.hartwig.hmftools.esvee.output.WriteType.READS;
import static com.hartwig.hmftools.esvee.output.WriteType.VCF;

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
import com.hartwig.hmftools.esvee.types.Junction;
import com.hartwig.hmftools.esvee.output.WriteType;
import com.hartwig.hmftools.esvee.read.Read;
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

    public final String VcfFile;
    public final List<WriteType> WriteTypes;

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
    public final boolean SkipDiscordant;
    public final boolean LogPhaseGroupLinks;
    public final int PhaseProcessingLimit;

    public final int Threads;

    public final String TruthsetFile;

    public static final String OUTPUT_VCF = "output_vcf";
    public static final String REF_GENOME_IMAGE = "ref_genome_image";
    public static final String DECOY_GENOME = "decoy_genome";
    public static final String JUNCTION_FILES = "junction_files";
    public static final String SPECIFIC_JUNCTIONS = "specific_junctions";

    public static final String WRITE_TYPES = "write_types";
    public static final String PERF_LOG_TIME = "perf_log_time";

    public static final String SKIP_DISCORDANT = "skip_discordant";
    public static final String PHASE_PROCESSING_LIMIT = "phase_process_limit";
    public static final String LOG_PHASE_GROUP_LINKS = "phase_group_links";
    public static final String RUN_LOCAL_PHASING = "run_local_phasing";

    public static final String ASSEMBLY_REF_BASE_WRITE_MAX = "asm_ref_base_write_max";

    public static final Logger SV_LOGGER = LogManager.getLogger(AssemblyConfig.class);

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

        JunctionFiles = Lists.newArrayList();

        if(configBuilder.hasValue(JUNCTION_FILES))
        {
            Arrays.stream(configBuilder.getValue(JUNCTION_FILES).split(CONFIG_FILE_DELIM)).forEach(x -> JunctionFiles.add(x));
        }
        else
        {
            String bamPath = pathFromFile(TumorBams.get(0));
            List<String> combinedSampleIds = Lists.newArrayList(TumorIds);
            combinedSampleIds.addAll(ReferenceIds);

            for(String sampleId : combinedSampleIds)
            {
                String junctionFile = bamPath + sampleId + SV_PREP_JUNCTIONS_FILE_ID;

                if(Files.exists(Paths.get(junctionFile)))
                    JunctionFiles.add(junctionFile);
            }
        }

        WriteTypes = Lists.newArrayList();

        if(configBuilder.hasValue(WRITE_TYPES))
        {
            String writeTypesStr = configBuilder.getValue(WRITE_TYPES);

            if(writeTypesStr.equals(WriteType.ALL))
            {
                Arrays.stream(WriteType.values()).filter(x -> x != READS).forEach(x -> WriteTypes.add(x));
            }
            else
            {
                String[] writeTypes = writeTypesStr.split(ITEM_DELIM, -1);
                Arrays.stream(writeTypes).forEach(x -> WriteTypes.add(WriteType.valueOf(x)));
            }
        }
        else
        {
            WriteTypes.add(VCF);
            WriteTypes.add(ASSEMBLY_BAM);
        }

        BamToolPath = configBuilder.getValue(BAMTOOL_PATH);

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        RefGenome = loadRefGenome(RefGenomeFile);

        RefGenomeImageFile = configBuilder.hasValue(REF_GENOME_IMAGE) ?
                configBuilder.getValue(REF_GENOME_IMAGE) : RefGenomeFile + REF_GENOME_IMAGE_EXTENSION;

        DecoyGenome = configBuilder.getValue(DECOY_GENOME);

        BamStringency = ValidationStringency.STRICT;

        RefGenomeCoords = RefGenVersion == V37 ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;

        if(!configBuilder.hasValue(OUTPUT_VCF) && !configBuilder.hasValue(OUTPUT_DIR))
        {
            SV_LOGGER.error("VCF output file or output directory required config");
            System.exit(1);
        }

        String vcfFile = configBuilder.getValue(OUTPUT_VCF);
        OutputDir = configBuilder.hasValue(OUTPUT_DIR) ? parseOutputDir(configBuilder) : pathFromFile(vcfFile);
        OutputId = configBuilder.getValue(OUTPUT_ID);

        VcfFile = vcfFile != null ? vcfFile : outputFilename(VCF);

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
        }

        if(WriteTypes.contains(READS) && !SpecificChrRegions.hasFilters() && SpecificJunctions.isEmpty())
        {
            SV_LOGGER.warn("writing assembly reads to TSV without region filtering may result in large output files & impact performance");
        }

        mLogReadIds = parseLogReadIds(configBuilder);
        mCheckLogReadIds = !mLogReadIds.isEmpty();

        SkipDiscordant = configBuilder.hasFlag(SKIP_DISCORDANT);
        PerfDebug = configBuilder.hasFlag(PERF_DEBUG);
        PerfLogTime = configBuilder.getDecimal(PERF_LOG_TIME);

        PhaseProcessingLimit = configBuilder.getInteger(PHASE_PROCESSING_LIMIT);

        AssemblyRefBaseWriteMax = configBuilder.getInteger(ASSEMBLY_REF_BASE_WRITE_MAX);
        LogPhaseGroupLinks = configBuilder.hasFlag(LOG_PHASE_GROUP_LINKS);

        Threads = parseThreads(configBuilder);

        TruthsetFile = TruthsetAnnotation.filename(configBuilder);
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

    public String outputFilename(final WriteType writeType)
    {
        return formOutputFile(OutputDir, TumorIds.get(0), writeType.fileId(), OutputId);
    }

    public void logReadId(final SAMRecord record, final String caller)
    {
        if(mCheckLogReadIds)
            logReadId(record.getReadName(), caller);
    }

    public void logReadId(final Read read, final String caller)
    {
        if(mCheckLogReadIds)
            logReadId(read.getName(), caller);
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

        configBuilder.addConfigItem(OUTPUT_VCF, false, "Output VCF filename");

        configBuilder.addPaths(
                JUNCTION_FILES, false, "List of SvPrep junction files, separated by ',', default is to match by sample name");

        addRefGenomeConfig(configBuilder, true);
        configBuilder.addPath(REF_GENOME_IMAGE, false, REFERENCE_BAM_DESC);
        configBuilder.addPath(DECOY_GENOME, false, "Decoy genome image file");

        String writeTypes = Arrays.stream(WriteType.values()).map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM));
        configBuilder.addConfigItem(WRITE_TYPES, false, "Write types from list: " + writeTypes);

        configBuilder.addConfigItem(LOG_READ_IDS, false, LOG_READ_IDS_DESC);
        configBuilder.addConfigItem(
                SPECIFIC_JUNCTIONS, false,
                "Specific junctions: format: Chromosome:Position:Orientation:Type (I, D or none), separated by ';'");

        configBuilder.addFlag(PERF_DEBUG, PERF_DEBUG_DESC);
        configBuilder.addDecimal(PERF_LOG_TIME, "Log performance data for routine exceeding specified time (0 = disabled)", 0);
        configBuilder.addFlag(SKIP_DISCORDANT, "Skip processing discordant-only groups");
        configBuilder.addFlag(LOG_PHASE_GROUP_LINKS, "Log assembly links to build phase groups");

        configBuilder.addInteger(
                ASSEMBLY_REF_BASE_WRITE_MAX, "Cap assembly ref bases in TSV and VCF, use zero to write all",
                DEFAULT_ASSEMBLY_REF_BASE_WRITE_MAX);

        configBuilder.addInteger(
                PHASE_PROCESSING_LIMIT, "Exclude phase groups above this size from extension and phase sets", 0);

        TruthsetAnnotation.registerConfig(configBuilder);
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

        BamStringency = ValidationStringency.SILENT;

        VcfFile = null;
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
        SkipDiscordant = true;
        LogPhaseGroupLinks = false;
        PhaseProcessingLimit = 0;
        Threads = 0;
        TruthsetFile = null;
    }
}
