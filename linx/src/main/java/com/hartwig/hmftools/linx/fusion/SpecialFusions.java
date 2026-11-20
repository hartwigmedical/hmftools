package com.hartwig.hmftools.linx.fusion;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.ENHANCER_KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.ENHANCER_PROMISCUOUS;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.NONE;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_ENHANCER_TARGET;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.ENHANCER;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_5P;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.ENHANCER_PROMISCUOUS_MIN_DISTANCE;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.ENHANCER_TRANSCRIPT_PREFIX;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.MAX_ENHANCER_PARTNER_GENE_DISTANCE_NO_ORIENT;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.MAX_UPSTREAM_DISTANCE_ENHANCER_KNOWN;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.MIN_ENHANCER_PARTNER_GENE_DISTANCE;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
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
        for(KnownFusionData kfData : mKnownFusionCache.getData())
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

                    for(ChrBaseRegion altRegion : kfData.getThreeGeneAltRegions())
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

        for(SvVarData var : svList)
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

        if(requiredOrientation != null && breakend.orientation() != requiredOrientation)
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

    public void checkEnhancerFusion(
            final BreakendGeneData startGene, final BreakendGeneData endGene, final List<GeneFusion> potentialFusions)
    {
        // These are allowed to fuse with 5’UTR splice donors from 1MB upstream. Phasing is assumed to be -1
        // In some specific cases, allow fusions in the 3’UTR region or even downstream of the gene (eg BCL2)
        // in which case the 3' gene's orientation must be facing back upstream into the coding region

        // orientation in the IG region is no longer checked
        boolean startIsEnhancerGene = mKnownFusionCache.withinEnhancerRegion(startGene.chromosome(), startGene.position());
        boolean endIsEnhancerGene = mKnownFusionCache.withinEnhancerRegion(endGene.chromosome(), endGene.position());

        if(startIsEnhancerGene == endIsEnhancerGene) // only one end can be
            return;

        BreakendGeneData enhancerGene = startIsEnhancerGene ? startGene : endGene;
        BreakendGeneData downGene = startIsEnhancerGene ? endGene : startGene;

        KnownFusionType knownType = NONE;

        KnownFusionData kfData = mKnownFusionCache.getDataByType(ENHANCER_KNOWN_PAIR).stream()
                .filter(x -> x.ThreeGene.equals(downGene.geneName()))
                .filter(x -> x.withinGeneRegion(enhancerGene.chromosome(), enhancerGene.position()))
                .findFirst().orElse(null);

        if(kfData == null)
        {
            kfData = mKnownFusionCache.getDataByType(ENHANCER_PROMISCUOUS).stream()
                    .filter(x -> x.withinGeneRegion(enhancerGene.chromosome(), enhancerGene.position()))
                    .findFirst().orElse(null);

            // check within the promiscuous region bounds
            if(kfData == null)
                return;

            if(enhancerGene.chromosome().equals(downGene.chromosome()))
            {
                int distance = abs(enhancerGene.position() - downGene.position());

                if(distance < MIN_ENHANCER_PARTNER_GENE_DISTANCE)
                    return;
            }
        }

        knownType = kfData.Type;

        List<BreakendTransData> downGeneTranscripts = Lists.newArrayList();

        int maxDownStreamDistance = -1;
        boolean allowCoding = false;

        if(kfData.geneStrand() == null)
        {
            maxDownStreamDistance = MAX_ENHANCER_PARTNER_GENE_DISTANCE_NO_ORIENT;
            allowCoding = true;
        }
        else
        {
            if(kfData.downstreamDistance(FS_DOWN) > 0)
                maxDownStreamDistance = kfData.downstreamDistance(FS_DOWN);
        }

        for(BreakendTransData transData : downGene.transcripts())
        {
            // typically must be prior to the 1st coding base of the 3’ gene partner
            if(transData.CodingBases == 0 && !transData.isUpstream())
            {
                if(transData.getDistanceUpstream() <= MAX_UPSTREAM_DISTANCE_ENHANCER_KNOWN)
                    downGeneTranscripts.add(transData);
            }
            else if(transData.codingType() == CODING && allowCoding)
            {
                downGeneTranscripts.add(transData);
            }
            else if(transData.codingType() == UTR_3P && transData.isUpstream() && maxDownStreamDistance > 0)
            {
                int downstreamDistance = transData.getDistanceDownstream();

                if(downstreamDistance <= maxDownStreamDistance)
                    downGeneTranscripts.add(transData);
            }
        }

        BreakendTransData upTrans = generateEnhancerTranscript(enhancerGene, kfData);

        // add new candidates only if no matching IG fusion has already been created
        // this is done since multiple upstream genes in the IG region can be processed, but all creating identical IG fusions
        for(BreakendTransData downTrans : downGeneTranscripts)
        {
            if(hasMatchingEnhancerFusion(potentialFusions, knownType, upTrans, downTrans))
                continue;

            GeneFusion fusion = new GeneFusion(upTrans, downTrans, true);
            fusion.setId(mFusionFinder.nextFusionId());
            fusion.setKnownType(knownType);
            potentialFusions.add(fusion);
        }
    }

    private static boolean hasMatchingEnhancerFusion(
            final List<GeneFusion> potentialFusions, final KnownFusionType knownType,
            final BreakendTransData upTrans, final BreakendTransData downTrans)
    {
        for(GeneFusion fusion : potentialFusions)
        {
            if(fusion.knownType() != knownType)
                continue;

            if(fusion.svId(true) != upTrans.breakendGeneData().varId() || fusion.svId(false) != downTrans.breakendGeneData().varId())
                continue;

            if(!fusion.upstreamTrans().TransData.TransName.equals(upTrans.TransData.TransName))
                continue;

            if(fusion.downstreamTrans() == downTrans) // straight comparison since not internally generated per fusion
                return true;
        }

        return false;
    }

    private BreakendTransData generateEnhancerTranscript(final BreakendGeneData gene, final KnownFusionData knownFusionData)
    {
        TranscriptData transData = new TranscriptData(
                0, String.format("%s%s", ENHANCER_TRANSCRIPT_PREFIX, knownFusionData.FiveGene), gene.geneId(), false,
                gene.strand(), knownFusionData.geneRegion().start(), knownFusionData.geneRegion().end(),
                null, null, "", null);

        BreakendTransData transcript = new BreakendTransData(
                gene, transData,  -1, -1, PHASE_NONE, PHASE_NONE, 0, 0);

        transcript.setCodingType(ENHANCER);
        transcript.setRegionType(TranscriptRegionType.IG);
        transcript.setIsDisruptive(true);
        return transcript;
    }

}
