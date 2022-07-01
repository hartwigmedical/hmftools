package com.hartwig.hmftools.isofox.common;

import static com.hartwig.hmftools.common.utils.ConfigUtils.loadGeneIdsFile;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_REGIONS;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_REGIONS_DESC;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.loadSpecificChromsomesOrRegions;
import static com.hartwig.hmftools.isofox.IsofoxConfig.GENE_ID_FILE;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.ENRICHED_GENE_BUFFER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.EXCLUDED_REGION_1_REF_37;
import static com.hartwig.hmftools.isofox.IsofoxConstants.EXCLUDED_REGION_1_REF_38;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ITEM_DELIM;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.SUB_ITEM_DELIM;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.isofox.IsofoxConstants;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class GeneRegionFilters
{
    public final List<String> SpecificChromosomes;
    public final List<ChrBaseRegion> SpecificRegions;

    public final List<String> RestrictedGeneIds; // specific set of genes to process
    public final List<String> EnrichedGeneIds; // genes to count by not fully process for any functional purpose
    public final ChrBaseRegion ExcludedRegion;

    public final List<ChrBaseRegion> RestrictedGeneRegions; // limit analysis to these regions only
    public final List<ChrBaseRegion> ImmuneGeneRegions;

    private final RefGenomeVersion mRefGenomeVersion;
    private final List<ChrBaseRegion> mExcludedGeneRegions; // exclude these regions based on geneId, enriched or excluded regions

    // config
    public static final String RESTRICTED_GENE_IDS = "restricted_gene_ids";
    private static final String ENRICHED_GENE_IDS = "enriched_gene_ids";

    public GeneRegionFilters(final RefGenomeVersion refGenomeVersion)
    {
        RestrictedGeneIds = Lists.newArrayList();
        EnrichedGeneIds = Lists.newArrayList();
        SpecificChromosomes = Lists.newArrayList();
        SpecificRegions = Lists.newArrayList();

        mExcludedGeneRegions = Lists.newArrayList();
        RestrictedGeneRegions = Lists.newArrayList();
        ImmuneGeneRegions = Lists.newArrayList();

        mRefGenomeVersion = refGenomeVersion;
        ExcludedRegion = refGenomeVersion.is37() ? EXCLUDED_REGION_1_REF_37 : EXCLUDED_REGION_1_REF_38;
    }

    public static void addCommandLineOptions(final Options options)
    {
        options.addOption(RESTRICTED_GENE_IDS, true, "Optional list of Ensmebl GeneIds separated by ';'");
        options.addOption(ENRICHED_GENE_IDS, true, "Optional list of geneIds to treat as enriched");
        addSpecificChromosomesRegionsConfig(options);
        options.addOption(SPECIFIC_REGIONS, true, SPECIFIC_REGIONS_DESC);
    }

    public void loadConfig(final CommandLine cmd) throws Exception
    {
        if(cmd.hasOption(ENRICHED_GENE_IDS))
        {
            Arrays.stream(cmd.getOptionValue(ENRICHED_GENE_IDS).split(ITEM_DELIM)).forEach(x -> EnrichedGeneIds.add(x));
        }
        else
        {
            IsofoxConstants.populateEnrichedGeneIds(EnrichedGeneIds, mRefGenomeVersion);
        }

        IsofoxConstants.populateImmuneRegions(ImmuneGeneRegions, mRefGenomeVersion);

        if(cmd.hasOption(GENE_ID_FILE))
        {
            final String inputFile = cmd.getOptionValue(GENE_ID_FILE);
            RestrictedGeneIds.addAll(loadGeneIdsFile(inputFile));

            if(!RestrictedGeneIds.isEmpty())
            {
                ISF_LOGGER.info("file({}) loaded {} restricted genes", inputFile, RestrictedGeneIds.size());
            }
        }
        else if(cmd.hasOption(RESTRICTED_GENE_IDS))
        {
            RestrictedGeneIds.addAll(Arrays.stream(cmd.getOptionValue(RESTRICTED_GENE_IDS).split(ITEM_DELIM)).collect(Collectors.toList()));
        }

        try
        {
            loadSpecificChromsomesOrRegions(cmd, SpecificChromosomes, SpecificRegions, ISF_LOGGER);
        }
        catch(ParseException e)
        {
        }
    }

    public boolean excludeChromosome(final String chromosome)
    {
        return !SpecificChromosomes.isEmpty() && !SpecificChromosomes.contains(chromosome);
    }

    public boolean skipRead(final String chromosome, int position)
    {
        // currently only used to filter out chimeric reads
        if(!HumanChromosome.contains(chromosome))
            return true;

        if(!SpecificChromosomes.isEmpty() && !SpecificChromosomes.contains(chromosome))
            return true;

        if(!SpecificRegions.isEmpty() && !SpecificRegions.stream().anyMatch(x -> x.containsPosition(chromosome, position)))
            return true;

        if(!RestrictedGeneRegions.isEmpty() && !RestrictedGeneRegions.stream().anyMatch(x -> x.containsPosition(chromosome, position)))
            return true;

        if(mExcludedGeneRegions.stream().anyMatch(x -> x.containsPosition(chromosome, position)))
            return true;

        return false;
    }

    public void buildGeneRegions(final EnsemblDataCache geneTransCache)
    {
        mExcludedGeneRegions.add(ExcludedRegion);

        EnrichedGeneIds.stream()
                .map(x -> geneTransCache.getGeneDataById(x))
                .filter(x -> x != null)
                .forEach(x -> mExcludedGeneRegions.add(new ChrBaseRegion(
                        x.Chromosome, x.GeneStart - ENRICHED_GENE_BUFFER, x.GeneEnd + ENRICHED_GENE_BUFFER)));

        RestrictedGeneIds.stream()
                .map(x -> geneTransCache.getGeneDataById(x))
                .filter(x -> x != null)
                .forEach(x -> RestrictedGeneRegions.add(new ChrBaseRegion(
                        x.Chromosome, x.GeneStart - 1000, x.GeneEnd + 1000)));
    }
}
