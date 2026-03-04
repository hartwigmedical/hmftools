package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.ENHANCER_KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_ENHANCER_TARGET;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.ENHANCER;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_5P;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.ENHANCER_PROMISCUOUS_MIN_DISTANCE;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.MAX_UPSTREAM_DISTANCE_ENHANCER_KNOWN;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.gene.TranscriptRegionType;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.linx.gene.BreakendGeneData;
import com.hartwig.hmftools.linx.gene.BreakendTransData;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvVarData;

public class SpecialFusions
{
    private final EnsemblDataCache mGeneDataCache;
    private final FusionFinder mFusionFinder;
    private final KnownFusionCache mKnownFusionCache;

    public SpecialFusions(final EnsemblDataCache ensemblDataCache, final FusionFinder fusionFinder)
    {
        mGeneDataCache = ensemblDataCache;
        mFusionFinder = fusionFinder;
        mKnownFusionCache = fusionFinder.getKnownFusionCache();
    }

    public void clear() {}

    public Map<String,Integer> getSpecificPreGeneDistances()
    {
        Map<String,Integer> specificGeneDistances = Maps.newHashMap();
        List<KnownFusionData> enhancerKnownPairs = mKnownFusionCache.getDataByType(ENHANCER_KNOWN_PAIR);
        enhancerKnownPairs.forEach(x -> specificGeneDistances.put(x.ThreeGene, MAX_UPSTREAM_DISTANCE_ENHANCER_KNOWN));
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

        List<KnownFusionData> enhancerTargets = mKnownFusionCache.getDataByType(PROMISCUOUS_ENHANCER_TARGET);

        if(enhancerTargets.isEmpty())
            return fusions;

        for(final SvVarData var : svList)
        {
            if(var.isSglBreakend())
                continue;

            for(int se = SE_START; se <= SE_END; ++se)
            {
                SvBreakend breakend = var.getBreakend(se);
                final int seIndex = se;

                enhancerTargets.stream()
                        .map(x -> formEnhancerTargetFusion(x, var, breakend, seIndex))
                        .filter(x -> x != null)
                        .forEach(x -> fusions.add(x));
            }

        }

        return fusions;
    }

    private GeneFusion formEnhancerTargetFusion(
            final KnownFusionData knownFusionData, final SvVarData var, final SvBreakend breakend, int seIndex)
    {
        Orientation requiredOrientation = knownFusionData.geneStrand();

        if(requiredOrientation != null && breakend.orientation() != requiredOrientation.asByte())
            return null;

        if(!knownFusionData.withinGeneRegion(breakend.chromosome(), breakend.position()))
            return null;

        if(var.type() != BND && var.length() < ENHANCER_PROMISCUOUS_MIN_DISTANCE)
            return null;

        SvBreakend otherBreakend = var.getBreakend(switchIndex(seIndex));

        // use any upstream breakend if it exists
        BreakendTransData upTrans = null;
        boolean otherIsStart = seIndex == SE_END;
        List<BreakendGeneData> upstreamGenes = var.getGenesList(otherIsStart);

        upTrans = upstreamGenes.stream().map(x -> x.canonical()).filter(x -> x != null).findFirst().orElse(null);

        if(upTrans == null)
        {
            GeneData geneData = new GeneData(
                    "", "", otherBreakend.chromosome(), POS_STRAND,
                    otherBreakend.position(), otherBreakend.position(), "");

            BreakendGeneData gene = new BreakendGeneData(var.id(), otherIsStart, geneData);
            gene.setSvData(var.getSvData(), var.jcn());

            TranscriptData transData = new TranscriptData(
                    0, "", "", false, POS_STRAND, 0, 0, null, null, "", null);

            BreakendTransData transcript = new BreakendTransData(
                    gene, transData,  0, 0, PHASE_NONE, PHASE_NONE, 0, 0);

            transcript.setCodingType(ENHANCER);
            transcript.setRegionType(TranscriptRegionType.UNKNOWN);

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

            GeneData geneData = mFusionFinder.getGeneTransCache().getGeneDataByName(knownFusionData.ThreeGene);
            TranscriptData transData = null;

            if(geneData != null)
            {
                transData = mFusionFinder.getGeneTransCache().getCanonicalTranscriptData(geneData.GeneId);
            }
            else
            {
                // note that these are configured to have the 'gene strand' being the required orientation of the 3' breakend
                // whereas the gene strand is actually the opposite
                byte geneStrand = (knownFusionData.geneStrand() == null || knownFusionData.geneStrand().isForward()) ? NEG_STRAND : POS_STRAND;

                geneData = new GeneData(
                        "", knownFusionData.ThreeGene, geneRegion.chromosome(), geneStrand,
                        geneRegion.start(), geneRegion.end(), "");

                transData = new TranscriptData(
                        0, "", knownFusionData.ThreeGene, false, geneStrand,
                        geneRegion.start(), geneRegion.end(), null, null, "", null);
            }

            BreakendGeneData gene = new BreakendGeneData(var.id(), seIndex == SE_START, geneData);
            gene.setSvData(var.getSvData(), var.jcn());

            BreakendTransData transcript = new BreakendTransData(
                    gene, transData,  1, 1, PHASE_NONE, PHASE_NONE, 0, 0);

            transcript.setCodingType(UTR_5P);
            transcript.setRegionType(TranscriptRegionType.UPSTREAM);

            downTrans = transcript;
        }

        GeneFusion fusion = new GeneFusion(upTrans, downTrans, false);
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
}
