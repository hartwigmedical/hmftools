package com.hartwig.hmftools.esvee;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
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
import static com.hartwig.hmftools.esvee.SvConstants.REF_GENOME_IMAGE_EXTENSION;
import static com.hartwig.hmftools.esvee.SvConstants.SV_PREP_JUNCTIONS_FILE_ID;
import static com.hartwig.hmftools.esvee.WriteType.ASSEMBLY_BAM;
import static com.hartwig.hmftools.esvee.WriteType.ASSEMBLY_READS;
import static com.hartwig.hmftools.esvee.WriteType.VCF;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.read.Read;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;

public class SvConfig
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

    public final ValidationStringency BamStringency;

    public final String VcfFile;
    public final List<WriteType> WriteTypes;

    public final String OutputDir;
    public final String OutputId;

    public final SpecificRegions SpecificChrRegions;

    public final boolean PerfDebug;
    public final double PerfLogTime;
    public final boolean OtherDebug;
    private final List<String> mLogReadIds;
    private final boolean mCheckLogReadIds;

    public final int Threads;

    public static final String OUTPUT_VCF = "output_vcf";
    public static final String REF_GENOME_IMAGE = "ref_genome_image";
    public static final String JUNCTION_FILES = "junction_files";

    // alternative to specifically load a tumor and/or ref sample and BAM
    public static final String SAMPLE_IDS = "samples";
    public static final String TUMOR_IDS = "samples";
    public static final String REFERENCE_IDS = "samples";
    public static final String SAMPLE_BAMS = "bam_files";
    public static final String TUMOR_BAMS = "bam_files";
    public static final String REFERENCE_BAMS = "bam_files";

    public static final String WRITE_TYPES = "write_types";
    public static final String HTML_SUMMARY_DIR = "html_dir";
    public static final String PLOT_DIAGRAMS = "plot_diagrams";
    public static final String OTHER_DEBUG = "other_debug";
    public static final String PERF_LOG_TIME = "perf_log_time";

    public static final Logger SV_LOGGER = LogManager.getLogger(SvConfig.class);

    public SvConfig(final ConfigBuilder configBuilder)
    {
        TumorIds = Arrays.stream(configBuilder.getValue(TUMOR).split(CONFIG_FILE_DELIM)).collect(Collectors.toList());
        ReferenceIds = Arrays.stream(configBuilder.getValue(REFERENCE).split(CONFIG_FILE_DELIM)).collect(Collectors.toList());

        TumorBams = Arrays.stream(configBuilder.getValue(TUMOR_BAM).split(CONFIG_FILE_DELIM)).collect(Collectors.toList());
        ReferenceBams = Arrays.stream(configBuilder.getValue(REFERENCE_BAM).split(CONFIG_FILE_DELIM)).collect(Collectors.toList());

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
                Arrays.stream(WriteType.values()).filter(x -> x != ASSEMBLY_READS).forEach(x -> WriteTypes.add(x));
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

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        RefGenome = loadRefGenome(RefGenomeFile);

        RefGenomeImageFile = configBuilder.hasValue(REF_GENOME_IMAGE) ?
                configBuilder.getValue(REF_GENOME_IMAGE) : RefGenomeFile + REF_GENOME_IMAGE_EXTENSION;

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

        mLogReadIds = parseLogReadIds(configBuilder);
        mCheckLogReadIds = !mLogReadIds.isEmpty();

        PerfDebug = configBuilder.hasFlag(PERF_DEBUG);
        PerfLogTime = configBuilder.getDecimal(PERF_LOG_TIME);
        OtherDebug = configBuilder.hasFlag(OTHER_DEBUG);

        Threads = parseThreads(configBuilder);
    }

    public String tumorBam() { return TumorBams.get(0); }
    public String referenceBam() { return !ReferenceBams.isEmpty() ? ReferenceBams.get(1) : null; }

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
        String filename = OutputDir;

        filename += TumorIds.get(0);

        if(OutputId != null)
            filename += "." + OutputId;

        filename += "." + writeType.fileId();

        return filename;
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
        configBuilder.addConfigItem(HTML_SUMMARY_DIR, false, "Directory for HTML summaries, default 'html'");
        configBuilder.addFlag(PLOT_DIAGRAMS, "Create HTML files containing SVGs");

        String writeTypes = Arrays.stream(WriteType.values()).map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM));
        configBuilder.addConfigItem(WRITE_TYPES, false, "Write types from list: " + writeTypes);

        configBuilder.addConfigItem(LOG_READ_IDS, false, LOG_READ_IDS_DESC);
        configBuilder.addFlag(PERF_DEBUG, PERF_DEBUG_DESC);
        configBuilder.addDecimal(PERF_LOG_TIME, "Log performance data for routine exceeding specified time (0 = disabled)", 0);
        configBuilder.addFlag(OTHER_DEBUG, "Various other debugging");

        SpecificRegions.addSpecificChromosomesRegionsConfig(configBuilder);

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
    }
}
