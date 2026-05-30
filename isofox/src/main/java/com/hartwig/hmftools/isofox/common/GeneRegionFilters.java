package com.hartwig.hmftools.isofox.common;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.loadChrBaseRegions;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.GENE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.GENE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadGeneIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.ENRICHED_GENE_BUFFER;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.isofox.IsofoxConstants;

import htsjdk.samtools.SAMRecord;

public class GeneRegionFilters
{
    public SpecificRegions SpecificChrRegions;

    public final List<String> RestrictedGeneIds; // limit expression analysis to a set of panel genes
    public final List<String> EnrichedGeneIds; // genes to count by not fully process for any functional purpose
    public final Map<String,List<BaseRegion>> ExcludedRegions;

    public final List<ChrBaseRegion> ImmuneGeneRegions;

    private final RefGenomeVersion mRefGenomeVersion;
    private boolean mHasSpecificRegions;

    // config
    private static final String ENRICHED_GENE_IDS = "enriched_gene_ids";
    private static final String EXCLUDED_REGIONS = "excluded_regions";

    public GeneRegionFilters(final RefGenomeVersion refGenomeVersion)
    {
        RestrictedGeneIds = Lists.newArrayList();
        EnrichedGeneIds = Lists.newArrayList();
        SpecificChrRegions = new SpecificRegions();
        mHasSpecificRegions = false;

        ImmuneGeneRegions = Lists.newArrayList();
        ExcludedRegions = Maps.newHashMap();

        mRefGenomeVersion = refGenomeVersion;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(GENE_ID_FILE, false, GENE_ID_FILE_DESC);
        configBuilder.addConfigItem(ENRICHED_GENE_IDS, "List of geneIds to treat as enriched");
        configBuilder.addPath(EXCLUDED_REGIONS, false, "List of excluded regions");
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

        if(configBuilder.hasValue(EXCLUDED_REGIONS))
        {
            String excludedRegionsFile = configBuilder.getValue(EXCLUDED_REGIONS);
            ExcludedRegions.putAll(loadChrBaseRegions(excludedRegionsFile, false));
            ISF_LOGGER.info("file({}) loaded {} excluded regions", excludedRegionsFile, ExcludedRegions.size());
        }

        SpecificChrRegions = SpecificRegions.from(configBuilder);
        mHasSpecificRegions = SpecificChrRegions.hasFilters();
    }

    public boolean excludeChromosome(final String chromosome)
    {
        return SpecificChrRegions.excludeChromosome(chromosome);
    }

    private static final int READ_END_BUFFER = 150;

    public boolean skipRead(final SAMRecord read, boolean checkMateAndSupp)
    {
        if(skipRead(read.getContig(), read.getAlignmentStart(), read.getAlignmentEnd()))
            return true;

        if(checkMateAndSupp)
        {
            // simple, non-cigar aware read end
            if(!read.getMateUnmappedFlag())
            {
                int mateReadStart = read.getMateAlignmentStart();
                if(skipRead(read.getMateReferenceName(), mateReadStart, mateReadStart))
                {
                    return true;
                }
            }

            SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(read);
            if(suppData != null)
            {
                if(skipRead(suppData.Chromosome, suppData.Position, suppData.Position))
                    return true;
            }
        }

        return false;
    }

    public boolean skipRead(final String chromosome, int readStart)
    {
        return skipRead(chromosome, readStart, readStart + READ_END_BUFFER);
    }

    public boolean skipRead(final String chromosome, int readStart, int readEnd)
    {
        // currently only used to filter out chimeric reads
        if(!HumanChromosome.contains(chromosome))
            return true;

        if(mHasSpecificRegions)
        {
            if(SpecificChrRegions.excludeChromosome(chromosome))
                return true;

            if(SpecificChrRegions.Regions.stream().noneMatch(x -> x.overlaps(chromosome, readStart, readEnd)))
                return true;
        }

        List<BaseRegion> excludedRegions = ExcludedRegions.get(chromosome);

        if(excludedRegions != null)
        {
            if(excludedRegions.stream().anyMatch(x -> positionsOverlap(x.start(), x.end(), readStart, readEnd)))
                return true;
        }

        return false;
    }

    public List<BaseRegion> findExcludedRegions(final ChrBaseRegion region)
    {
        List<BaseRegion> excludedRegions = ExcludedRegions.get(region.Chromosome);

        if(excludedRegions == null)
            return Collections.emptyList();

        return excludedRegions.stream().filter(x -> x.overlaps(region)).collect(Collectors.toList());
    }

    public void buildGeneRegions(final EnsemblDataCache geneTransCache)
    {
        // add the regions from any enriched genes to the excluded regions set
        for(String enrichedGeneId : EnrichedGeneIds)
        {
            GeneData geneData = geneTransCache.getGeneDataById(enrichedGeneId);

            if(geneData == null)
            {
                ISF_LOGGER.warn("enriched gene ID({}) missing from Ensembl cache", enrichedGeneId);
                continue;
            }

            List<BaseRegion> excludedRegions = ExcludedRegions.get(geneData.Chromosome);

            if(excludedRegions == null)
            {
                excludedRegions = Lists.newArrayList();
                ExcludedRegions.put(geneData.Chromosome, excludedRegions);
            }

            excludedRegions.add(new BaseRegion(
                    geneData.GeneStart - ENRICHED_GENE_BUFFER, geneData.GeneEnd + ENRICHED_GENE_BUFFER));
        }

        for(List<BaseRegion> excludedRegions : ExcludedRegions.values())
        {
            Collections.sort(excludedRegions);
        }
    }
}
