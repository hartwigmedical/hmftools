package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.StructuralVariant;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class GermlineVcfConfig
{
    // run config
    public final String SampleId;
    public final String OutputDir;
    public final String Scope;
    public final String VcfFile;
    public final String OutputVcfFile;
    public final String VcfsFile;
    public final String ProcessedFile;
    public final String BatchRunRootDir;
    public final boolean LinkByAssembly;

    // filtering config
    public final boolean RequireGridssPass;
    public final boolean LogFiltered;
    public final int QualScoreThreshold;
    public final List<String> RestrictedChromosomes;
    public final boolean RequireGene;

    private static final String SAMPLE = "sample";
    private static final String VCF_FILE = "vcf";
    private static final String VCFS_FILE = "vcfs_file";
    private static final String BATCH_ROOT_DIR = "batch_root_dir";
    private static final String OUTPUT_VCF = "output_vcf";

    private static final String PROCESSED_FILE = "processed";

    private static final String SCOPE = "scope";
    public static final String GENE_PANEL_FILE = "gene_panel_file";
    private static final String LINK_BY_ASSEMBLY = "link_by_assembly";
    private static final String CHECK_DISRUPTIONS = "check_disruptions";

    private static final String REQUIRE_PASS = "require_pass";
    private static final String LOG_FILTERED = "log_filtered";
    private static final String QUAL_SCORE_THRESHOLD = "qs_threshold";
    private static final String SPECIFIC_CHROMOSOMES = "specific_chr";
    private static final String REQUIRE_GENE = "require_gene";

    public GermlineVcfConfig(final CommandLine cmd)
    {
        SampleId = cmd.getOptionValue(SAMPLE);
        OutputDir = parseOutputDir(cmd);

        VcfFile = cmd.getOptionValue(VCF_FILE, "");
        VcfsFile = cmd.getOptionValue(VCFS_FILE, "");
        BatchRunRootDir = cmd.getOptionValue(BATCH_ROOT_DIR, "");
        OutputVcfFile = cmd.getOptionValue(VCF_FILE, "");

        RestrictedChromosomes = cmd.hasOption(SPECIFIC_CHROMOSOMES) ?
                Arrays.stream(cmd.getOptionValue(SPECIFIC_CHROMOSOMES, "")
                .split(";")).collect(Collectors.toList()) : Lists.newArrayList();

        // unused
        RequireGridssPass = cmd.hasOption(REQUIRE_PASS);
        QualScoreThreshold = Integer.parseInt(cmd.getOptionValue(QUAL_SCORE_THRESHOLD, "350"));
        ProcessedFile = cmd.getOptionValue(PROCESSED_FILE, "");
        Scope = cmd.getOptionValue(SCOPE);
        LinkByAssembly = cmd.hasOption(LINK_BY_ASSEMBLY);

        LogFiltered = cmd.hasOption(LOG_FILTERED);
        RequireGene = cmd.hasOption(REQUIRE_GENE);
    }

    public boolean excludeVariant(final StructuralVariant sv)
    {
        // optionally filter out all but specified chromosomes
        if(!RestrictedChromosomes.isEmpty() && !RestrictedChromosomes.contains(sv.chromosome(true))
        && (sv.type() == SGL || !RestrictedChromosomes.contains(sv.chromosome(false))))
        {
            return true;
        }

        if(RequireGridssPass && !sv.filter().contains(PASS))
        {
            return true;
        }


        return false;
    }

    public static void addCommandLineOptions(final Options options)
    {
        options.addOption(SAMPLE, true, "Name of the tumor sample");
        options.addOption(VCF_FILE, true, "Path to the GRIDSS structural variant VCF file");
        options.addOption(BATCH_ROOT_DIR, true, "Path to the root directory for sample runs");
        options.addOption(PROCESSED_FILE, true, "Path to a previously-run output file");
        options.addOption(ENSEMBL_DATA_DIR, true, "Ensembl data cache directory");
        options.addOption(GENE_PANEL_FILE, true, "Gene panel file");
        options.addOption(CHECK_DISRUPTIONS, false, "Check gene disruptions and filter out non-disruptive genes");
        options.addOption(LINK_BY_ASSEMBLY, false, "Look for assembled links");
        options.addOption(SCOPE, true, "Scope: germline or somatic");
        options.addOption(OUTPUT_DIR, true, "Path to write results");
        options.addOption(OUTPUT_VCF, true, "Path to write results");
        options.addOption(LOG_DEBUG, false, "Log verbose");
        options.addOption(REF_GENOME, true, REF_GENOME_CFG_DESC);

        options.addOption(REQUIRE_PASS, false, "Require variants to have GRIDSS filter = PASS");
        options.addOption(QUAL_SCORE_THRESHOLD, true, "Qual score threshold");
        options.addOption(SPECIFIC_CHROMOSOMES, true, "Optional set of chromosomes to restrict search to");
        options.addOption(LOG_FILTERED, false, "Log filtered variants");
        options.addOption(REQUIRE_GENE, false, "Only log SVs linked to a gene panel entry");

        PonCache.addCmdLineArgs(options);
        HotspotCache.addCmdLineArgs(options);
    }
}
