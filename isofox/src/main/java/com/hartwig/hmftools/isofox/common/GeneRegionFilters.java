package com.hartwig.hmftools.isofox.common;

import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.GENE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.GENE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadGeneIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
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

    public final List<ChrBaseRegion> ImmuneGeneRegions;

    private final RefGenomeVersion mRefGenomeVersion;
    private boolean mHasSpecificRegions;

    // config
    private static final String ENRICHED_GENE_IDS = "enriched_gene_ids";

    public GeneRegionFilters(final RefGenomeVersion refGenomeVersion)
    {
        RestrictedGeneIds = Lists.newArrayList();
        EnrichedGeneIds = Lists.newArrayList();
        SpecificChrRegions = new SpecificRegions();
        mHasSpecificRegions = false;

        ImmuneGeneRegions = Lists.newArrayList();

        mRefGenomeVersion = refGenomeVersion;
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

        SpecificChrRegions = SpecificRegions.from(configBuilder);
        mHasSpecificRegions = SpecificChrRegions.hasFilters();
    }

    public boolean excludeChromosome(final String chromosome)
    {
        return SpecificChrRegions.excludeChromosome(chromosome);
    }

    private static final int READ_END_BUFFER = 150; // NOTE: could use BAM sampling read length

    public boolean skipRead(final SAMRecord read, boolean checkMateAndSupp)
    {
        if(skipRead(read.getContig(), read.getAlignmentStart(), read.getAlignmentEnd(), false))
            return true;

        if(checkMateAndSupp)
        {
            // simple, non-cigar aware read end
            if(!read.getMateUnmappedFlag())
            {
                int mateReadStart = read.getMateAlignmentStart();
                if(skipRead(read.getMateReferenceName(), mateReadStart, mateReadStart, true))
                {
                    return true;
                }
            }

            SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(read);
            if(suppData != null)
            {
                if(skipRead(suppData.Chromosome, suppData.Position, suppData.Position, true))
                    return true;
            }
        }

        return false;
    }

    public boolean skipRead(final String chromosome, int readStart)
    {
        return skipRead(chromosome, readStart, readStart + READ_END_BUFFER, false);
    }

    public boolean skipRead(final String chromosome, int readStart, int readEnd, boolean isMateOrSupp)
    {
        // currently only used to filter out chimeric reads
        if(!HumanChromosome.contains(chromosome))
            return true;

        if(!isMateOrSupp && mHasSpecificRegions)
        {
            if(SpecificChrRegions.excludeChromosome(chromosome))
                return true;

            if(SpecificChrRegions.Regions.stream().noneMatch(x -> x.overlaps(chromosome, readStart, readEnd)))
                return true;
        }

        return false;
    }
}
