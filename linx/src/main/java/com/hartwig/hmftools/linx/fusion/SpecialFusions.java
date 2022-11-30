package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_ENHANCER_TARGET;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.ENHANCER;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_5P;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.MAX_UPSTREAM_DISTANCE_IG_KNOWN;
import static com.hartwig.hmftools.linx.fusion.FusionReportability.determineReportability;
import static com.hartwig.hmftools.linx.fusion.ReportableReason.OK;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.gene.TranscriptRegionType;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.linx.gene.BreakendGeneData;
import com.hartwig.hmftools.linx.gene.BreakendTransData;
import com.hartwig.hmftools.linx.types.SglMapping;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

public class SpecialFusions
{
    private final EnsemblDataCache mGeneDataCache;
    private final FusionFinder mFusionFinder;
    private final FusionConfig mFusionConfig;
    private final KnownFusionCache mKnownFusionCache;

    public SpecialFusions(final EnsemblDataCache ensemblDataCache, final FusionFinder fusionFinder, final FusionConfig fusionConfig)
    {
        mGeneDataCache = ensemblDataCache;
        mFusionFinder = fusionFinder;
        mKnownFusionCache = fusionFinder.getKnownFusionCache();
        mFusionConfig = fusionConfig;
    }

    public void clear() {}

    public Map<String,Integer> getSpecificPreGeneDistances()
    {
        Map<String,Integer> specificGeneDistances = Maps.newHashMap();
        List<KnownFusionData> igKnownPairs = mKnownFusionCache.getDataByType(IG_KNOWN_PAIR);
        igKnownPairs.forEach(x -> specificGeneDistances.put(x.ThreeGene, MAX_UPSTREAM_DISTANCE_IG_KNOWN));
        return specificGeneDistances;
    }

    public void cacheSpecialFusionGenes()
    {
        for(final KnownFusionData kfData : mKnownFusionCache.getData())
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

    public List<GeneFusion> findFusions(final List<SvVarData> svList)
    {
        List<GeneFusion> fusions = Lists.newArrayList();

        // always report SVs by themselves
        final List<SvVarData> knownPairSglCandidates = Lists.newArrayList();

        List<KnownFusionData> enhancerTargets = mKnownFusionCache.getDataByType(PROMISCUOUS_ENHANCER_TARGET);

        if(enhancerTargets.isEmpty())
            return fusions;

        for(final SvVarData var : svList)
        {
            // if(var.isSglBreakend() && !var.getGenesList(true).isEmpty())
            //    knownPairSglCandidates.add(var);

            for(int se = SE_START; se <= SE_END; ++se)
            {
                SvBreakend breakend = var.getBreakend(se);
                final int seIndex = se;

                if(var.isSglBreakend() && se == SE_END)
                {
                    for(SglMapping mapping : var.getSglMappings())
                    {
                        enhancerTargets.stream()
                                .map(x -> formEnhancerTargetFusion(x, var, null, mapping, seIndex))
                                .filter(x -> x != null)
                                .forEach(x -> fusions.add(x));
                    }
                }
                else
                {
                    enhancerTargets.stream()
                            .map(x -> formEnhancerTargetFusion(x, var, breakend, null, seIndex))
                            .filter(x -> x != null)
                            .forEach(x -> fusions.add(x));
                }
            }

        }

        /*
        if(!knownPairSglCandidates.isEmpty())
            fusions.addAll(findKnownPairSglFusions(knownPairSglCandidates));
        */

        return fusions;
    }

    private GeneFusion formEnhancerTargetFusion(
            final KnownFusionData knownFusionData, final SvVarData var, final SvBreakend breakend, final SglMapping sglMapping, int seIndex)
    {
        if(sglMapping != null)
        {
            if(!knownFusionData.matchesGeneRegion(sglMapping.Chromosome, sglMapping.Position, sglMapping.Orientation))
                return null;

        }
        else
        {
            if(!knownFusionData.matchesGeneRegion(breakend.chromosome(), breakend.position(), breakend.orientation()))
                return null;
        }

        final SvBreakend otherBreakend = var.getBreakend(switchIndex(seIndex));

        // use any upstream breakend if it exists
        BreakendTransData upTrans = null;
        boolean otherIsStart = seIndex == SE_END;
        List<BreakendGeneData> upstreamGenes = var.getGenesList(otherIsStart);

        if(!upstreamGenes.isEmpty())
        {
            BreakendGeneData gene = upstreamGenes.get(0);
            upTrans = gene.canonical();
        }
        else
        {
            String otherChr;
            int otherPosition;

            if(otherBreakend != null)
            {
                otherChr = otherBreakend.chromosome();
                otherPosition = otherBreakend.position();
            }
            else
            {
                otherChr = !var.getSglMappings().isEmpty() ? var.getSglMappings().get(0).Chromosome : "";
                otherPosition = !var.getSglMappings().isEmpty() ? var.getSglMappings().get(0).Position : 0;
            }

            GeneData geneData = new GeneData(
                    "", "", otherChr, POS_STRAND, otherPosition, otherPosition, "");

            BreakendGeneData gene = new BreakendGeneData(var.id(), otherIsStart, geneData);
            gene.setSvData(var.getSvData(), var.jcn());

            TranscriptData transData = new TranscriptData(
                    0, "", "", false, POS_STRAND, 0, 0, null, null, "");

            BreakendTransData transcript = new BreakendTransData(
                    gene, transData,  0, 0, PHASE_NONE, PHASE_NONE, 0, 0);

            transcript.setCodingType(ENHANCER);
            transcript.setRegionType(TranscriptRegionType.UNKNOWN);
            // transcript.setIsDisruptive(true);

            upTrans = transcript;
        }

        BreakendTransData downTrans = null;
        BreakendGeneData downstreamGene = var.getGenesList(seIndex == SE_START).stream()
                .filter(x -> x.geneName().equals(knownFusionData.ThreeGene))
                .findFirst().orElse(null);

        if(downstreamGene != null)
        {
            downTrans = downstreamGene.canonical();
        }
        else
        {
            ChrBaseRegion geneRegion = knownFusionData.geneRegion();
            GeneData geneData = new GeneData(
                    "", knownFusionData.ThreeGene, geneRegion.chromosome(), knownFusionData.geneStrand(),
                    geneRegion.start(), geneRegion.end(), "");

            BreakendGeneData gene = new BreakendGeneData(var.id(), seIndex == SE_START, geneData);
            gene.setSvData(var.getSvData(), var.jcn());

            TranscriptData transData = new TranscriptData(
                    0, "", knownFusionData.ThreeGene, false, POS_STRAND, geneRegion.start(), geneRegion.end(),
                    null, null, "");

            BreakendTransData transcript = new BreakendTransData(
                    gene, transData,  1, 1, PHASE_NONE, PHASE_NONE, 0, 0);

            transcript.setCodingType(UTR_5P);
            transcript.setRegionType(TranscriptRegionType.UPSTREAM);
            // transcript.setIsDisruptive(true);

            downTrans = transcript;
        }

        GeneFusion fusion = new GeneFusion(upTrans, downTrans, false);
        // fusion.setReportable(true);
        fusion.setId(mFusionFinder.nextFusionId());
        fusion.setKnownType(PROMISCUOUS_ENHANCER_TARGET);

        FusionAnnotations annotations = ImmutableFusionAnnotations.builder()
                .clusterId(var.getCluster().id())
                .clusterCount(var.getCluster().getSvCount())
                .resolvedType(var.getCluster().getResolvedType().toString())
                .chainInfo(null)
                .terminatedUp(false)
                .terminatedDown(false)
                .build();

        fusion.setAnnotations(annotations);

        return fusion;
    }

    private List<GeneFusion> findKnownPairSglFusions(final List<SvVarData> knownPairSglCandidates)
    {
        List<GeneFusion> fusions = Lists.newArrayList();

        if(knownPairSglCandidates.size() < 2)
            return fusions;

        List<KnownFusionData> knownPairData = mKnownFusionCache.knownPairData();

        for(int i = 0; i < knownPairSglCandidates.size() - 1; ++i)
        {
            SvVarData sgl1 = knownPairSglCandidates.get(i);

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

            for(int j = i + 1; j < knownPairSglCandidates.size(); ++j)
            {
                SvVarData sgl2 = knownPairSglCandidates.get(j);

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
