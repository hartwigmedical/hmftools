package com.hartwig.hmftools.isofox.common;

import static com.hartwig.hmftools.common.region.ExcludedRegions.getPolyGRegion;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.loadSpecificChromsomesOrRegions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.GENE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.GENE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadGeneIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.ENRICHED_GENE_BUFFER;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.isofox.IsofoxConstants;

import org.apache.commons.cli.ParseException;

import htsjdk.samtools.SAMRecord;

public class GeneRegionFilters
{
    public final SpecificRegions SpecificChrRegions;

    public final List<String> RestrictedGeneIds; // specific set of genes to process
    public final List<String> EnrichedGeneIds; // genes to count by not fully process for any functional purpose
    public final ChrBaseRegion ExcludedRegion;

    public final List<ChrBaseRegion> RestrictedGeneRegions; // limit analysis to these regions only
    public final List<ChrBaseRegion> ImmuneGeneRegions;

    private final RefGenomeVersion mRefGenomeVersion;
    private final List<ChrBaseRegion> mExcludedGeneRegions; // exclude these regions based on geneId, enriched or excluded regions

    // config
    private static final String ENRICHED_GENE_IDS = "enriched_gene_ids";

    public GeneRegionFilters(final RefGenomeVersion refGenomeVersion)
    {
        RestrictedGeneIds = Lists.newArrayList();
        EnrichedGeneIds = Lists.newArrayList();
        SpecificChrRegions = new SpecificRegions();

        mExcludedGeneRegions = Lists.newArrayList();
        RestrictedGeneRegions = Lists.newArrayList();
        ImmuneGeneRegions = Lists.newArrayList();

        mRefGenomeVersion = refGenomeVersion;
        ExcludedRegion = getPolyGRegion(refGenomeVersion);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(GENE_ID_FILE, false, GENE_ID_FILE_DESC);
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

        try
        {
            loadSpecificChromsomesOrRegions(configBuilder, SpecificChrRegions.Chromosomes, SpecificChrRegions.Regions);
        }
        catch(ParseException e)
        {
            System.exit(1);
        }
    }

    public boolean excludeChromosome(final String chromosome)
    {
        return SpecificChrRegions.excludeChromosome(chromosome);
    }

    public boolean skipRead(final String chromosome, int position)
    {
        // currently only used to filter out chimeric reads
        if(!HumanChromosome.contains(chromosome))
            return true;

        if(SpecificChrRegions.excludeChromosome(chromosome))
            return true;

        if(SpecificChrRegions.excludePosition(chromosome, position))
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
