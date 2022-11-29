package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_KNOWN_PAIR;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.MAX_UPSTREAM_DISTANCE_IG_KNOWN;
import static com.hartwig.hmftools.linx.fusion.FusionReportability.determineReportability;
import static com.hartwig.hmftools.linx.fusion.ReportableReason.OK;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.linx.gene.BreakendGeneData;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

public class SpecialFusions
{
    private final EnsemblDataCache mGeneDataCache;
    private final FusionFinder mFusionFinder;
    private final FusionConfig mFusionConfig;
    private final List<SvVarData> mKnownPairSglCandidates;

    public SpecialFusions(final EnsemblDataCache ensemblDataCache, final FusionFinder fusionFinder, final FusionConfig fusionConfig)
    {
        mGeneDataCache = ensemblDataCache;
        mFusionFinder = fusionFinder;
        mFusionConfig = fusionConfig;
        mKnownPairSglCandidates = Lists.newArrayList();
    }

    public void clear()
    {
        mKnownPairSglCandidates.clear();
    }

    public Map<String,Integer> getSpecificPreGeneDistances()
    {
        Map<String,Integer> specificGeneDistances = Maps.newHashMap();
        List<KnownFusionData> igKnownPairs = mFusionFinder.getKnownFusionCache().getDataByType(IG_KNOWN_PAIR);
        igKnownPairs.forEach(x -> specificGeneDistances.put(x.ThreeGene, MAX_UPSTREAM_DISTANCE_IG_KNOWN));
        return specificGeneDistances;
    }

    public void cacheSpecialFusionGenes()
    {
        for(final KnownFusionData kfData : mFusionFinder.getKnownFusionCache().getData())
        {
            if(kfData.downstreamDistance(FS_UP) > 0)
            {
                final GeneData geneData = mGeneDataCache.getGeneDataByName(kfData.FiveGene);
                if(geneData != null)
                {
                    mGeneDataCache.addDownstreamGeneAnnotations(geneData, kfData.downstreamDistance(FS_UP));
                }
            }

            if(kfData.downstreamDistance(FS_DOWN) > 0)
            {
                final GeneData geneData = mGeneDataCache.getGeneDataByName(kfData.ThreeGene);
                if(geneData != null)
                {
                    mGeneDataCache.addDownstreamGeneAnnotations(geneData, kfData.downstreamDistance(FS_DOWN));
                }
            }

            if(!kfData.getThreeGeneAltRegions().isEmpty())
            {
                final GeneData geneData = mGeneDataCache.getGeneDataByName(kfData.ThreeGene);

                if(geneData != null)
                {
                    if(mGeneDataCache.getAlternativeGeneData().stream().anyMatch(x -> x.GeneId.equals(geneData.GeneId)))
                        continue;

                    for(final ChrBaseRegion altRegion : kfData.getThreeGeneAltRegions())
                    {
                        mGeneDataCache.getAlternativeGeneData().add(new GeneData(
                                geneData.GeneId, geneData.GeneName, altRegion.Chromosome, geneData.Strand, altRegion.start(), altRegion.end(), ""));
                    }
                }
            }
        }
    }

    public void findFusions(final List<SvVarData> svList)
    {
        // always report SVs by themselves
        for(final SvVarData var : svList)
        {
            if(var.isSglBreakend())
            {
                if(!var.getGenesList(true).isEmpty())
                    mKnownPairSglCandidates.add(var);

                if(var.getSglMappings().isEmpty())
                    continue;
            }
        }
    }

    public List<GeneFusion> findKnownPairSglFusions()
    {
        List<GeneFusion> fusions = Lists.newArrayList();

        if(mKnownPairSglCandidates.size() < 2)
            return fusions;

        List<KnownFusionData> knownPairData = mFusionFinder.getKnownFusionCache().knownPairData();

        for(int i = 0; i < mKnownPairSglCandidates.size() - 1; ++i)
        {
            SvVarData sgl1 = mKnownPairSglCandidates.get(i);

            final List<BreakendGeneData> genesList1 = sgl1.getGenesList(true);

            List<KnownFusionData> fivePrimeData = Lists.newArrayList();
            List<KnownFusionData> threePrimeData = Lists.newArrayList();

            for(BreakendGeneData geneData : genesList1)
            {
                if(geneData.isUpstream())
                {
                    knownPairData.stream()
                            .filter(x -> genesList1.stream().anyMatch(y -> y.geneName().equals(x.FiveGene)))
                            .forEach(x -> fivePrimeData.add(x));
                }
                else
                {
                    knownPairData.stream()
                            .filter(x -> genesList1.stream().anyMatch(y -> y.geneName().equals(x.ThreeGene)))
                            .forEach(x -> threePrimeData.add(x));
                }
            }

            for(int j = i + 1; j < mKnownPairSglCandidates.size(); ++j)
            {
                SvVarData sgl2 = mKnownPairSglCandidates.get(j);

                final List<BreakendGeneData> genesList2 = sgl2.getGenesList(true);

                for(BreakendGeneData geneData : genesList2)
                {
                    if(geneData.isUpstream())
                    {
                        for(KnownFusionData kpData : threePrimeData)
                        {
                            if(kpData.FiveGene.equals(geneData.geneName()))
                            {
                                checkKnownPairSgls(kpData, sgl1, sgl2);

                                fusions.addAll(mFusionFinder.findFusions(
                                        Lists.newArrayList(geneData),
                                        Lists.newArrayList(genesList1.stream()
                                                .filter(x -> x.geneName().equals(kpData.ThreeGene))
                                                .collect(Collectors.toList()))));
                            }
                        }
                    }
                    else
                    {
                        for(KnownFusionData kpData : fivePrimeData)
                        {
                            if(kpData.ThreeGene.equals(geneData.geneName()))
                            {
                                checkKnownPairSgls(kpData, sgl1, sgl2);

                                fusions.addAll(mFusionFinder.findFusions(
                                        Lists.newArrayList(genesList1.stream()
                                                .filter(x -> x.geneName().equals(kpData.FiveGene))
                                                .collect(Collectors.toList())),
                                        Lists.newArrayList(geneData)));
                            }
                        }
                    }
                }

                if(fusions.isEmpty())
                    continue;

                if(mFusionConfig.LogReportableOnly)
                {
                    fusions = fusions.stream().filter(x -> determineReportability(x) == OK).collect(Collectors.toList());
                }

                final SvCluster cluster = sgl1.getCluster();

                // check transcript disruptions
                for(final GeneFusion fusion : fusions)
                {
                    FusionAnnotations annotations = ImmutableFusionAnnotations.builder()
                            .clusterId(cluster.id())
                            .clusterCount(cluster.getSvCount())
                            .resolvedType(cluster.getResolvedType().toString())
                            .chainInfo(null)
                            .terminatedUp(false)
                            .terminatedDown(false)
                            .build();

                    fusion.setAnnotations(annotations);
                    fusions.add(fusion);
                }
            }
        }

        return fusions;
    }

    private boolean checkKnownPairSgls(final KnownFusionData kpData, final SvVarData sgl1, final SvVarData sgl2)
    {
        /*
        LNX_LOGGER.info("sample({}) sgls({} & {}) knownPair({}-{}) sequences({}, {})",
            mSampleId, sgl1.posId(), sgl2.posId(), kpData.FiveGene, kpData.ThreeGene,
            sgl1.getSvData().insertSequence(), sgl2.getSvData().insertSequence());
        */

        return false;
    }
}
