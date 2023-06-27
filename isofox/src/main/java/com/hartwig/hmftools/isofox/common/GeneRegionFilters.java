package com.hartwig.hmftools.isofox.common;

import static com.hartwig.hmftools.common.sv.ExcludedRegions.getPolyGRegion;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadGeneIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.loadSpecificChromsomesOrRegions;
import static com.hartwig.hmftools.isofox.IsofoxConfig.GENE_ID_FILE;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.ENRICHED_GENE_BUFFER;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.isofox.IsofoxConstants;

import org.apache.commons.cli.ParseException;

import htsjdk.samtools.SAMRecord;

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
        ExcludedRegion = getPolyGRegion(refGenomeVersion);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(RESTRICTED_GENE_IDS, "List of Ensmebl GeneIds separated by ';'");
        configBuilder.addConfigItem(ENRICHED_GENE_IDS, "List of geneIds to treat as enriched");
        addSpecificChromosomesRegionsConfig(configBuilder);
    }

    public void loadConfig(final ConfigBuilder configBuilder)
    {
        if(configBuilder.hasValue(ENRICHED_GENE_IDS))
        {
            Arrays.stream(configBuilder.getValue(ENRICHED_GENE_IDS).split(ITEM_DELIM)).forEach(x -> EnrichedGeneIds.add(x));
        }
        else
        {
            IsofoxConstants.populateEnrichedGeneIds(EnrichedGeneIds, mRefGenomeVersion);
        }

        IsofoxConstants.populateImmuneRegions(ImmuneGeneRegions, mRefGenomeVersion);

        if(configBuilder.hasValue(GENE_ID_FILE))
        {
            final String inputFile = configBuilder.getValue(GENE_ID_FILE);
            RestrictedGeneIds.addAll(loadGeneIdsFile(inputFile));

            if(!RestrictedGeneIds.isEmpty())
            {
                ISF_LOGGER.info("file({}) loaded {} restricted genes", inputFile, RestrictedGeneIds.size());
            }
        }
        else if(configBuilder.hasValue(RESTRICTED_GENE_IDS))
        {
            RestrictedGeneIds.addAll(Arrays.stream(configBuilder.getValue(RESTRICTED_GENE_IDS).split(ITEM_DELIM)).collect(Collectors.toList()));
        }

        try
        {
            loadSpecificChromsomesOrRegions(configBuilder, SpecificChromosomes, SpecificRegions, ISF_LOGGER);
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

    public static boolean inExcludedRegion(final ChrBaseRegion excludedRegion, final SAMRecord record)
    {
        if(excludedRegion == null)
            return false;

        if(excludedRegion.containsPosition(record.getStart()) || excludedRegion.containsPosition(record.getEnd()))
            return true;

        // check the mate as well
        if(record.getMateReferenceName() != null && excludedRegion.containsPosition(record.getMateReferenceName(), record.getMateAlignmentStart()))
        {
            return true;
        }

        return false;
    }
}
