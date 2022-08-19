package com.hartwig.hmftools.svprep.depth;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.loadSpecificRegions;
import static com.hartwig.hmftools.svprep.SvCommon.DELIM;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class DepthConfig
{
    public final String InputVcf;
    public final String OutputVcf;
    public final List<String> Samples;
    public final List<String> BamFiles;
    public final String RefGenome;
    public final RefGenomeVersion RefGenVersion;
    public final double VafCap;

    public final int Threads;
    public final List<ChrBaseRegion> SpecificRegions;

    private static final String INPUT_VCF = "input_vcf";
    private static final String OUTPUT_VCF = "output_vcf";
    private static final String SAMPLES = "samples";
    private static final String BAM_FILES = "bam_files";
    private static final String THREADS = "threads";
    private static final String VAF_CAP = "vaf_cap";

    protected static final int MAX_GAP = 500;
    protected static final double DEFAULT_VAF_CAP = 0.001;

    public DepthConfig(final CommandLine cmd)
    {
        InputVcf = cmd.getOptionValue(INPUT_VCF);
        OutputVcf = cmd.getOptionValue(OUTPUT_VCF);

        Samples = Arrays.stream(cmd.getOptionValue(SAMPLES).split(DELIM, -1)).collect(Collectors.toList());
        BamFiles = Arrays.stream(cmd.getOptionValue(BAM_FILES).split(DELIM, -1)).collect(Collectors.toList());

        RefGenome = cmd.getOptionValue(REF_GENOME);
        RefGenVersion = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION, V37.toString()));;
        VafCap = Double.parseDouble(cmd.getOptionValue(VAF_CAP, String.valueOf(DEFAULT_VAF_CAP)));
        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "1"));

        SpecificRegions = Lists.newArrayList();

        try
        {
            SpecificRegions.addAll(loadSpecificRegions(cmd));
        }
        catch(ParseException e)
        {
            SV_LOGGER.error("failed to load specific regions");
        }
    }

    public static void addOptions(final Options options)
    {
        options.addOption(INPUT_VCF, true, "Input VCF File");
        options.addOption(SAMPLES, true, "Sample IDs corresponding to BAM files");
        options.addOption(BAM_FILES, true, "BAM file(s) to slice for depth");
        options.addOption(OUTPUT_VCF, true, "Output VCF File");
        options.addOption(REF_GENOME, true, REF_GENOME_CFG_DESC);
        options.addOption(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);
        options.addOption(THREADS, true, "Multi-thread count");
        addSpecificChromosomesRegionsConfig(options);
    }
}
