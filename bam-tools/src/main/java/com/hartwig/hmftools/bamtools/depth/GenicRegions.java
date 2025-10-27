package com.hartwig.hmftools.bamtools.depth;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.depth.HighDepthCombiner.PANEL_HIGH_DEPTH_THRESHOLD;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.KnownFusionCache.KNOWN_FUSIONS_FILE;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.immune.ImmuneRegions.getIgRegion;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class GenicRegions
{
    private final RefGenomeVersion mRefGenVersion;
    private final List<DriverGene> mDriverGenes;
    private final KnownFusionCache mKnownFusionCache;
    private final EnsemblDataCache mEnsemblDataCache;
    private final boolean mRemoveGeneOverlaps;
    private final Map<String,List<BaseRegion>> mFixedGeneRegions;

    protected static final String FIXED_GENE_REGIONS = "fixed_gene_regions"; // panel genic regions, move any overlaps with these
    protected static final String REMOVE_GENE_OVERLAPS = "remove_gene_overlaps";

    public GenicRegions(final ConfigBuilder configBuilder)
    {
        mRefGenVersion = RefGenomeVersion.from(configBuilder);
        mRemoveGeneOverlaps = configBuilder.hasFlag(REMOVE_GENE_OVERLAPS);

        mKnownFusionCache = new KnownFusionCache();
        mKnownFusionCache.loadFromFile(configBuilder.getValue(KNOWN_FUSIONS_FILE));

        mDriverGenes = DriverGenePanelConfig.loadDriverGenes(configBuilder);

        if(configBuilder.hasValue(ENSEMBL_DATA_DIR))
        {
            mEnsemblDataCache = new EnsemblDataCache(configBuilder);
            mEnsemblDataCache.load(true);
        }
        else
        {
            mEnsemblDataCache = null;
        }

        mFixedGeneRegions = Maps.newHashMap();

        if(configBuilder.hasValue(FIXED_GENE_REGIONS))
        {
            mFixedGeneRegions.putAll(ChrBaseRegion.loadChrBaseRegions(configBuilder.getValue(FIXED_GENE_REGIONS)));

            if(!mFixedGeneRegions.isEmpty())
            {
                BT_LOGGER.debug("loaded {} fixed gene regions", mFixedGeneRegions.values().stream().mapToInt(x -> x.size()).sum());
            }
        }

        if(!mFixedGeneRegions.isEmpty() || !mDriverGenes.isEmpty())
        {
            BT_LOGGER.debug("OVERLAP,GeneName,Chromosome,GeneStart,GeneEnd,RegionStart,RegionEnd,OverlapBases,OverlapType,SampleCount,DepthMin,DepthMax");
        }
    }

    public void checkFixedGenicRegionOverlaps(final Map<String,List<HighDepthRegion>> finalRegions)
    {
        if(mFixedGeneRegions.isEmpty())
            return;

        for(Map.Entry<String,List<HighDepthRegion>> entry : finalRegions.entrySet())
        {
            String chromosome = entry.getKey();
            List<BaseRegion> geneRegions = mFixedGeneRegions.get(entry.getKey());

            if(geneRegions == null || geneRegions.isEmpty())
                continue;

            List<HighDepthRegion> highDepthRegions = entry.getValue();

            int index = 0;

            while(index < highDepthRegions.size())
            {
                HighDepthRegion highDepthRegion = highDepthRegions.get(index);

                BaseRegion genicRegion = geneRegions.stream().filter(x -> x.overlaps(highDepthRegion)).findFirst().orElse(null);

                if(genicRegion != null)
                {
                    int overlappingBases = min(genicRegion.end(), highDepthRegion.end()) - max(genicRegion.start(), highDepthRegion.start());

                    BT_LOGGER.trace("OVERLAP,{},{},{},{},{},{},{},{},{},{},{}",
                            chromosome, chromosome, genicRegion.start(), genicRegion.end(),
                            highDepthRegion.start(), highDepthRegion.end(), overlappingBases, "FIXED_GENE",
                            highDepthRegion.SampleCount, highDepthRegion.DepthMin, highDepthRegion.DepthMax);

                    highDepthRegions.remove(index);
                }
                else
                {
                    ++index;
                }
            }
        }
    }

    public void checkKnownGeneOverlaps(final Map<String,List<HighDepthRegion>> finalRegions)
    {
        if(mEnsemblDataCache == null || (mDriverGenes.isEmpty() && !mKnownFusionCache.hasValidData()))
            return;

        // convert driver and fusion genes into regions to compare
        Set<String> addedGenes = Sets.newHashSet();

        mDriverGenes.forEach(x -> checkGeneRegion(finalRegions, x.gene(), addedGenes, false));

        List<KnownFusionData> knownFusionData = mKnownFusionCache.getDataByType(KNOWN_PAIR);

        if(knownFusionData != null)
        {
            knownFusionData.forEach(x -> checkGeneRegion(finalRegions, x.FiveGene, addedGenes, false));
            knownFusionData.forEach(x -> checkGeneRegion(finalRegions, x.ThreeGene, addedGenes, true));
        }

        // check the IG regions
        checkGeneRegion(finalRegions, "IGH", addedGenes, false);
        checkGeneRegion(finalRegions, "IGK", addedGenes, false);
        checkGeneRegion(finalRegions, "IGL", addedGenes, false);
    }

    private void checkGeneRegion(
            final Map<String,List<HighDepthRegion>> finalRegions, final String geneName, final Set<String> addedGenes, boolean checkUpstream)
    {
        if(addedGenes.contains(geneName))
            return;

        addedGenes.add(geneName);

        GeneData geneData = mEnsemblDataCache.getGeneDataByName(geneName);
        if(geneData == null)
        {
            ChrBaseRegion igRegion = getIgRegion(geneName, mRefGenVersion);

            if(igRegion == null)
                return;

            geneData = new GeneData(geneName, geneName, igRegion.chromosome(), POS_STRAND, igRegion.start(), igRegion.end(), "");
        }

        List<HighDepthRegion> highDepthRegions = finalRegions.get(geneData.Chromosome);

        if(highDepthRegions == null)
            return;

        int geneStartPos = checkUpstream && geneData.Strand == POS_STRAND ?
                geneData.GeneStart - 10000 : geneData.GeneStart;

        int geneEndPos = checkUpstream && geneData.Strand == NEG_STRAND ? geneData.GeneEnd + 10000 : geneData.GeneEnd;

        int index = 0;
        while(index < highDepthRegions.size())
        {
            HighDepthRegion highDepthRegion = highDepthRegions.get(index);

            if(!positionsOverlap(geneStartPos, geneEndPos, highDepthRegion.start(), highDepthRegion.end()))
            {
                ++index;
                continue;
            }

            int overlappingBases = max(min(geneData.GeneEnd, highDepthRegion.end()) - max(geneData.GeneStart, highDepthRegion.start()) + 1, 0);

            String overlapType;

            if(overlappingBases == 0)
                overlapType = "UPSTREAM";
            else if(positionsWithin(highDepthRegion.start(), highDepthRegion.end(), geneStartPos, geneEndPos))
                overlapType = "WITHIN_GENE";
            else if(positionsWithin(geneStartPos, geneEndPos, highDepthRegion.start(), highDepthRegion.end()))
                overlapType = "CONTAINS_GENE";
            else
                overlapType = "OVERLAPS";

            BT_LOGGER.debug("OVERLAP,{},{},{},{},{},{},{},{},{},{},{}",
                    geneData.GeneName, geneData.Chromosome, geneData.GeneStart, geneData.GeneEnd,
                    highDepthRegion.start(), highDepthRegion.end(), overlappingBases, overlapType,
                    highDepthRegion.SampleCount, highDepthRegion.DepthMin, highDepthRegion.DepthMax);

            boolean removeGeneOverlaps = mRemoveGeneOverlaps && highDepthRegion.DepthMax < PANEL_HIGH_DEPTH_THRESHOLD;

            if(!removeGeneOverlaps)
            {
                ++index;
                continue;
            }

            // truncate or remove the region
            if(positionsWithin(highDepthRegion.start(), highDepthRegion.end(), geneStartPos, geneEndPos)
                    || positionsWithin(geneStartPos, geneEndPos, highDepthRegion.start(), highDepthRegion.end()))
            {
                highDepthRegions.remove(index);
                continue;
            }

            if(highDepthRegion.start() < geneStartPos)
                highDepthRegion.setEnd(geneStartPos - 1);
            else
                highDepthRegion.setStart(geneEndPos + 1);

            ++index;
        }
    }
}
