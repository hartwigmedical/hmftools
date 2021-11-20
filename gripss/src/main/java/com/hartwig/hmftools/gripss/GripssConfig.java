package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.sv.StructuralVariant;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class GripssConfig
{
    // run config
    public final String SampleId;
    public final String ReferenceId;
    public final RefGenomeVersion RefGenVersion;
    public final String VcfFile;

    public final String OutputDir;
    // public final String Scope;
    public final List<String> RestrictedChromosomes;

    private static final String SAMPLE = "sample";
    private static final String REFERENCE = "reference";
    private static final String VCF_FILE = "vcf";

    private static final String SCOPE = "scope";

    private static final String SPECIFIC_CHROMOSOMES = "specific_chr";

    public GripssConfig(final CommandLine cmd)
    {
        SampleId = cmd.getOptionValue(SAMPLE);
        ReferenceId = cmd.getOptionValue(REFERENCE);
        OutputDir = parseOutputDir(cmd);

        VcfFile = cmd.getOptionValue(VCF_FILE, "");
        RefGenVersion = cmd.hasOption(REF_GENOME_VERSION) ? RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION)) : V37;

        RestrictedChromosomes = cmd.hasOption(SPECIFIC_CHROMOSOMES) ?
                Arrays.stream(cmd.getOptionValue(SPECIFIC_CHROMOSOMES, "")
                .split(";")).collect(Collectors.toList()) : Lists.newArrayList();

        // Scope = cmd.getOptionValue(SCOPE);
    }

    public GripssConfig(
            final String sampleId, final String referenceId, final RefGenomeVersion refGenVersion, final String vcfFile)
    {
        SampleId = sampleId;
        ReferenceId = referenceId;
        RefGenVersion = refGenVersion;
        VcfFile = vcfFile;
        OutputDir = null;
        RestrictedChromosomes = Lists.newArrayList();
    }

    public boolean tumorOnly() { return ReferenceId.isEmpty(); }

    public boolean excludeVariant(final StructuralVariant sv)
    {
        // optionally filter out all but specified chromosomes
        if(!RestrictedChromosomes.isEmpty() && !RestrictedChromosomes.contains(sv.chromosome(true))
        && (sv.type() == SGL || !RestrictedChromosomes.contains(sv.chromosome(false))))
        {
            return true;
        }

        return false;
    }

    public static void addCommandLineOptions(final Options options)
    {
        options.addOption(SAMPLE, true, "Name of the tumor sample");
        options.addOption(REFERENCE, true, "Optional, name of the reference sample");
        options.addOption(VCF_FILE, true, "Path to the GRIDSS structural variant VCF file");
        options.addOption(REF_GENOME, true, REF_GENOME_CFG_DESC);
        options.addOption(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);
        options.addOption(SCOPE, true, "Scope: germline or somatic");

        addOutputDir(options);
        addLoggingOptions(options);
        addRefGenomeConfig(options);

        options.addOption(SPECIFIC_CHROMOSOMES, true, "Optional set of chromosomes to restrict search to");

        PonCache.addCmdLineArgs(options);
        HotspotCache.addCmdLineArgs(options);
    }
}
