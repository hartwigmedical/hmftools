package com.hartwig.hmftools.bamtools.slice;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.checkFileExists;
import static com.hartwig.hmftools.bamtools.metrics.MetricsConfig.BAM_FILE;
import static com.hartwig.hmftools.bamtools.metrics.MetricsConfig.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.metrics.MetricsConfig.LOG_READ_IDS;
import static com.hartwig.hmftools.bamtools.metrics.MetricsConfig.PARTITION_SIZE;
import static com.hartwig.hmftools.bamtools.metrics.MetricsConfig.PERF_DEBUG;
import static com.hartwig.hmftools.bamtools.metrics.MetricsConfig.SAMPLE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.DEFAULT_CHR_PARTITION_SIZE;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.loadSpecificChromsomesOrRegions;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bamtools.common.CommonUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class SliceConfig
{
    public final String SampleId;
    public final String BamFile;
    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;

    public final int PartitionSize;

    public final String OutputDir;
    public final String OutputId;

    public final boolean WriteBam;
    public final boolean WriteReads;
    public final int Threads;

    // debug
    public final List<String> SpecificChromosomes;
    public final List<ChrBaseRegion> SpecificRegions;
    public final boolean PerfDebug;

    private boolean mIsValid;

    private static final String SLICE_BED = "slice_bed";
    private static final String WRITE_BAM = "write_bam";
    private static final String WRITE_READS = "write_reads";

    public SliceConfig(final CommandLine cmd)
    {
        mIsValid = true;

        SampleId = cmd.getOptionValue(SAMPLE);
        BamFile = cmd.getOptionValue(BAM_FILE);
        RefGenomeFile = cmd.getOptionValue(REF_GENOME);
        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID);
        WriteReads = cmd.hasOption(WRITE_READS);
        WriteBam = cmd.hasOption(WRITE_BAM) || !WriteReads;

        if(SampleId == null || BamFile == null || OutputDir == null || RefGenomeFile == null)
        {
            BT_LOGGER.error("missing config: sample({}) bam({}) refGenome({}) outputDir({})",
                    SampleId != null, BamFile != null, RefGenomeFile != null, OutputDir != null);
            mIsValid = false;
        }

        RefGenVersion = cmd.hasOption(REF_GENOME_VERSION) ? RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION)) : V37;

        BT_LOGGER.info("refGenome({}), bam({})", RefGenVersion, BamFile);
        BT_LOGGER.info("output({})", OutputDir);

        PartitionSize = Integer.parseInt(cmd.getOptionValue(PARTITION_SIZE, String.valueOf(DEFAULT_CHR_PARTITION_SIZE)));

        SpecificChromosomes = Lists.newArrayList();
        SpecificRegions = Lists.newArrayList();

        if(cmd.hasOption(SLICE_BED))
        {
            loadSliceRegions(cmd.getOptionValue(SLICE_BED));
        }
        else
        {
            try
            {
                loadSpecificChromsomesOrRegions(cmd, SpecificChromosomes, SpecificRegions, BT_LOGGER);
                SpecificChromosomes.clear();
            }
            catch(ParseException e)
            {
                mIsValid = false;
            }
        }

        if(SpecificRegions.isEmpty())
        {
            BT_LOGGER.error("missing specific regions or slice BED file for slicing");
            mIsValid = false;
        }

        Threads = parseThreads(cmd);

        PerfDebug = cmd.hasOption(PERF_DEBUG);
    }

    public boolean isValid()
    {
        if(!mIsValid)
            return false;

        mIsValid = checkFileExists(BamFile) && checkFileExists(RefGenomeFile);
        return mIsValid;
    }

    public String formFilename(final String fileType) { return CommonUtils.formFilename(SampleId, OutputDir, OutputId, fileType); }

    private void loadSliceRegions(final String filename)
    {
        if(filename == null)
            return;

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));

            for(String line : lines)
            {
                if(line.contains("Chromosome"))
                    continue;

                final String[] values = line.split("\t", -1);

                if(values.length < 3)
                {
                    BT_LOGGER.error("invalid slice BED entry: {}", line);
                    mIsValid = false;
                    return;
                }

                String chromosome = values[0];
                int posStart = Integer.parseInt(values[1]) + 1;
                int posEnd = Integer.parseInt(values[2]);
                SpecificRegions.add(new ChrBaseRegion(chromosome, posStart, posEnd));
            }

            BT_LOGGER.info("loaded {} slice regions from file", SpecificRegions.size(), filename);
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to load hotspot data file({}): {}", filename, e.toString());
            mIsValid = false;
        }
    }

    public static Options createCmdLineOptions()
    {
        final Options options = new Options();
        addOutputOptions(options);
        addLoggingOptions(options);

        options.addOption(SAMPLE, true, "Tumor sample ID");
        options.addOption(BAM_FILE, true, "BAM file location");
        addRefGenomeConfig(options);;
        options.addOption(SLICE_BED, true, "BED file defining region to slice");
        options.addOption(PARTITION_SIZE, true, "Partition size, default: " + DEFAULT_CHR_PARTITION_SIZE);
        options.addOption(WRITE_BAM, false, "Write BAM file for sliced region");
        options.addOption(WRITE_READS, false, "Write CSV reads file for sliced region");
        addThreadOptions(options);

        addSpecificChromosomesRegionsConfig(options);
        options.addOption(LOG_READ_IDS, true, "Log specific read IDs, separated by ';'");
        options.addOption(PERF_DEBUG, false, "Detailed performance tracking and logging");

        return options;
    }
}
