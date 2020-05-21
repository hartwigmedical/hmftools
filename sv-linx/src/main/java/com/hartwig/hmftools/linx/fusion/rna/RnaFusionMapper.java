package com.hartwig.hmftools.linx.fusion.rna;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.EXON_RANK_MIN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWNSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UPSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.streamStr;
import static com.hartwig.hmftools.common.fusion.GeneFusion.REPORTABLE_TYPE_3P_PROM;
import static com.hartwig.hmftools.common.fusion.GeneFusion.REPORTABLE_TYPE_5P_PROM;
import static com.hartwig.hmftools.common.fusion.GeneFusion.REPORTABLE_TYPE_BOTH_PROM;
import static com.hartwig.hmftools.common.fusion.GeneFusion.REPORTABLE_TYPE_KNOWN;
import static com.hartwig.hmftools.common.fusion.GeneFusion.REPORTABLE_TYPE_NONE;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.linx.fusion.FusionFinder.checkFusionLogic;
import static com.hartwig.hmftools.common.fusion.KnownFusionData.FIVE_GENE;
import static com.hartwig.hmftools.common.fusion.KnownFusionData.THREE_GENE;
import static com.hartwig.hmftools.linx.fusion.rna.RnaFusionData.getRnaSourceDelimiter;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.fusion.GeneAnnotation;
import com.hartwig.hmftools.common.fusion.GeneFusion;
import com.hartwig.hmftools.common.fusion.Transcript;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.fusion.FusionFinder;
import com.hartwig.hmftools.linx.fusion.FusionParameters;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvVarData;

public class RnaFusionMapper
{
    private String mSampleId;
    private String mOutputDir;

    private final FusionFinder mFusionFinder;
    private final RnaMatchWriter mWriter;
    private final FusionParameters mFusionParams;
    private final EnsemblDataCache mGeneTransCache;
    private final Map<String, List<RnaFusionData>> mSampleRnaData;

    private final List<GeneFusion> mDnaFusions;
    private final Map<GeneFusion,String> mDnaInvalidFusions;

    public RnaFusionMapper(final String outputDir, final EnsemblDataCache geneTransCache, FusionFinder fusionFinder,
            final List<GeneFusion> dnaFusions, final Map<GeneFusion,String> dnaInvalidFusions)
    {
        mSampleRnaData = Maps.newHashMap();
        mFusionFinder = fusionFinder;
        mGeneTransCache = geneTransCache;
        mDnaFusions = dnaFusions;
        mDnaInvalidFusions = dnaInvalidFusions;
        mWriter = new RnaMatchWriter(outputDir);

        mFusionParams = new FusionParameters();
        mFusionParams.RequirePhaseMatch = false;
        mFusionParams.AllowExonSkipping = false;
    }

    public final Map<String, List<RnaFusionData>> getSampleRnaData() { return mSampleRnaData; }
    public final List<RnaFusionData> getSampleRnaData(final String sampleId) { return mSampleRnaData.get(sampleId); }

    public void assessRnaFusions(final String sampleId, Map<String, List<SvBreakend>> chrBreakendMap)
    {
        mSampleId = sampleId;

        final List<RnaFusionData> rnaFusionList = mSampleRnaData.get(mSampleId);

        if (rnaFusionList == null || rnaFusionList.isEmpty())
            return;

        LNX_LOGGER.debug("assessing {} RNA fusions", rnaFusionList.size());

        for (final RnaFusionData rnaFusion : rnaFusionList)
        {
            setRnaFusionData(rnaFusion);

            if(!rnaFusion.isValid())
                continue;

            annotateRnaFusions(rnaFusion, chrBreakendMap);

            mWriter.writeRnaMatchData(mSampleId, rnaFusion);
        }

        // move from consideration to de-link RNA data from SV types
        mSampleRnaData.remove(mSampleId);
    }

    private void annotateRnaFusions(final RnaFusionData rnaFusion, Map<String, List<SvBreakend>> chrBreakendMap)
    {
        /* Matching and annotation logic:
            - find all breakends in the RNA up and down gene
            - for them, find the any transcripts which a) have the exon boundary in the RNA position AND
            - b) are in the correct relative position:
                - upstream: at or after the RNA boundary down to the start of the next exon
                - downstream: at or before the RNA bounday up to the start of the preceding exon
                - If a transcript on the downstream gene starts on the 2nd exon, the fusion is allowed to match up to the nearer
                of a splice acceptor site with same orientation on a previous gene OR 100k bases upstream of the transcript.
                (is the distance up available for these??)
                - if multiple transcripts exist for the same breakend, take the canonical or longest (but it should make no difference)
            - if multiple breakends meet these criteria at either end, prioritise in the following order
                - both breakends are either end of the same structural variant
                - both breakends are in the same chain
                - both breakends are in the same cluster
                - otherwise take the nearest breakend to the RNA position
            - if no breakend is found on either upstream or downstream gene meeting the above criteria then record the nearest ID,
            distance and min number of skipped splice sites.
        */

        // viable breakends and their matching transcript
        List<SvBreakend> viableUpBreakends = Lists.newArrayList();
        List<SvBreakend> viableDownBreakends = Lists.newArrayList();
        List<Transcript> viableUpTranscripts = Lists.newArrayList();
        List<Transcript> viableDownTranscripts = Lists.newArrayList();

        // transcripts on the correct side and orientation of the RNA boundary
        List<Transcript> nearUpTranscripts = Lists.newArrayList();
        List<Transcript> nearDownTranscripts = Lists.newArrayList();
        List<SvBreakend> nearUpBreakends = Lists.newArrayList();
        List<SvBreakend> nearDownBreakends = Lists.newArrayList();

        // non-viable transcripts to be used if no others are found
        List<Transcript> genicUpTranscripts = Lists.newArrayList();
        List<Transcript> genicDownTranscripts = Lists.newArrayList();
        List<SvBreakend> genicUpBreakends = Lists.newArrayList();
        List<SvBreakend> genicDownBreakends = Lists.newArrayList();

        boolean isExactRnaExon = rnaFusion.matchesKnownSpliceSite();

        for(int fs = FS_UPSTREAM; fs <= FS_DOWNSTREAM ; ++fs)
        {
            boolean isUpstream = (fs == 0);
            final String chromosome = rnaFusion.Chromosomes[fs];
            int rnaPosition = rnaFusion.Positions[fs];
            byte geneStrand = rnaFusion.Strands[fs];
            List<SvBreakend> viableBreakends = isUpstream ? viableUpBreakends : viableDownBreakends;
            List<SvBreakend> nearBreakends = isUpstream ? nearUpBreakends : nearDownBreakends;
            List<SvBreakend> genicBreakends = isUpstream ? genicUpBreakends : genicDownBreakends;
            List<Transcript> viableTranscripts = isUpstream ? viableUpTranscripts : viableDownTranscripts;
            List<Transcript> nearTranscripts = isUpstream ? nearUpTranscripts : nearDownTranscripts;
            List<Transcript> genicTranscripts = isUpstream ? genicUpTranscripts : genicDownTranscripts;
            String geneName = rnaFusion.GeneNames[fs];

            final List<SvBreakend> breakendList = chrBreakendMap.get(chromosome);

            if(breakendList == null)
                continue;

            for(final SvBreakend breakend : breakendList)
            {
                final SvVarData var = breakend.getSV();

                if(var.isInferredSgl())
                    continue;

                // check whether breakend falls in genic region
                List<GeneAnnotation> genesList = var.getGenesList(breakend.usesStart())
                        .stream()
                        .filter(x -> x.GeneName.equals(geneName))
                        .collect(Collectors.toList());

                if(genesList.isEmpty())
                    continue;

                // check that breakend has correct orientation and position relative to RNA breakend
                boolean correctLocation = isViableBreakend(breakend, rnaPosition, geneStrand, isUpstream);

                // check whether any of the breakend's transcripts match the exon (exact or nearest) of the RNA fusion breakpoint
                for(final Transcript trans : genesList.get(0).transcripts())
                {
                    if(isExactRnaExon)
                    {
                        if(!rnaFusion.getExactMatchTransIds()[fs].contains(trans.StableId))
                            continue;
                    }
                    else if(!trans.isCanonical())
                    {
                        continue;
                    }

                    if(correctLocation)
                    {
                        nearBreakends.add(breakend);
                        nearTranscripts.add(trans);
                    }
                    else
                    {
                        genicBreakends.add(breakend);
                        genicTranscripts.add(trans);
                    }

                    if(correctLocation && isTranscriptBreakendViableForRnaBoundary(
                            trans, isUpstream,  breakend.position(), rnaPosition, isExactRnaExon))
                    {
                        viableBreakends.add(breakend);
                        viableTranscripts.add(trans);
                    }
                }
            }
        }

        LNX_LOGGER.debug("rna fusion({}) breakend matches: upstream(viable={} near={} genic={}) downstream(viable={} near={} genic={})",
                rnaFusion.name(), viableUpBreakends.size(), nearUpBreakends.size(), genicUpBreakends.size(),
                viableDownBreakends.size(), nearDownBreakends.size(), genicDownBreakends.size());

        // run them through fusion logic (ie a pair of breakend lists), but don't require phase matching
        if(!viableUpBreakends.isEmpty() && !viableDownBreakends.isEmpty())
        {
            GeneFusion topCandidateFusion = null;
            SvBreakend topUpBreakend = null;
            SvBreakend topDownBreakend = null;
            boolean topCandidateFusionViable = false;

            for (int i = 0; i < viableUpBreakends.size(); ++i)
            {
                final SvBreakend upBreakend = viableUpBreakends.get(i);
                final Transcript upTrans = viableUpTranscripts.get(i);

                if(upBreakend.getSV().isSglBreakend())
                    continue;

                for (int j = 0; j < viableDownBreakends.size(); ++j)
                {
                    final SvBreakend downBreakend = viableDownBreakends.get(j);
                    final Transcript downTrans = viableDownTranscripts.get(j);

                    if(downBreakend.getSV().isSglBreakend())
                        continue;

                    GeneFusion possibleFusion = checkFusionLogic(upTrans, downTrans, mFusionParams);
                    boolean viableFusion = possibleFusion != null;

                    // form one any way but mark it as not meeting standard fusion rules
                    if(possibleFusion == null)
                    {
                        possibleFusion = new GeneFusion(upTrans, downTrans, false);
                    }

                    if (topCandidateFusion == null
                    || isCandidateBetter(topCandidateFusion, topUpBreakend, topDownBreakend, possibleFusion, upBreakend, downBreakend,
                            rnaFusion, topCandidateFusionViable, viableFusion))
                    {
                        topCandidateFusion = possibleFusion;
                        topUpBreakend = upBreakend;
                        topDownBreakend = downBreakend;
                        topCandidateFusionViable = viableFusion;

                        LNX_LOGGER.debug("rnaFusion({}) first pair({} & {})", rnaFusion.name(), upBreakend.toString(), downBreakend.toString());
                    }
                }
            }

            if(topCandidateFusion != null)
            {
                rnaFusion.setTranscriptData(
                        FS_UPSTREAM, topCandidateFusion.upstreamTrans(), topUpBreakend,
                        true, true,  0);

                rnaFusion.setTranscriptData(
                        FS_DOWNSTREAM, topCandidateFusion.downstreamTrans(), topDownBreakend,
                        true, true,0);

                rnaFusion.setViableFusion(topCandidateFusionViable, topCandidateFusion.phaseMatched());
            }
        }
        else
        {
            // select the closest breakend's transcript
            for(int fs = FS_UPSTREAM; fs <= FS_DOWNSTREAM ; ++fs)
            {
                boolean isUpstream = (fs == 0);
                int rnaPosition = rnaFusion.Positions[fs];

                List<Transcript> transcriptList;
                List<SvBreakend> breakendList;
                boolean isViable = false;
                boolean correctLocation = false;

                // use the viable transcripts if present, otherwise the nearest
                if(isUpstream)
                {
                    if(!viableUpTranscripts.isEmpty())
                    {
                        isViable = true;
                        correctLocation = true;
                        transcriptList = viableUpTranscripts;
                        breakendList = viableUpBreakends;
                    }
                    else if(!nearUpTranscripts.isEmpty())
                    {
                        correctLocation = true;
                        transcriptList = nearUpTranscripts;
                        breakendList = nearUpBreakends;
                    }
                    else
                    {
                        transcriptList = genicUpTranscripts;
                        breakendList = genicUpBreakends;
                    }
                }
                else
                {
                    if(!viableDownTranscripts.isEmpty())
                    {
                        isViable = true;
                        correctLocation = true;
                        transcriptList = viableDownTranscripts;
                        breakendList = viableDownBreakends;
                    }
                    else if(!nearDownTranscripts.isEmpty())
                    {
                        correctLocation = true;
                        transcriptList = nearDownTranscripts;
                        breakendList = nearDownBreakends;
                    }
                    else
                    {
                        transcriptList = genicDownTranscripts;
                        breakendList = genicDownBreakends;
                    }
                }

                Transcript closestTrans = null;
                SvBreakend closestBreakend = null;
                int closestDistance = 0;

                for (int j = 0; j < transcriptList.size(); ++j)
                {
                    final Transcript trans = transcriptList.get(j);
                    final SvBreakend breakend = breakendList.get(j);

                    int distance = abs(rnaPosition - trans.svPosition());
                    if(closestTrans == null || distance < closestDistance)
                    {
                        closestDistance = distance;
                        closestTrans = trans;
                        closestBreakend = breakend;
                    }
                }

                if(closestTrans != null)
                {
                    int exonsSkipped = 0;

                    if(!isViable)
                    {
                        // for non-viable breakends, provide the exons skipped count
                        String geneId = closestTrans.gene().StableId;
                        final int rnaExonData[] = mGeneTransCache.getExonRankings(geneId, rnaPosition);
                        final int svPosExonData[] = mGeneTransCache.getExonRankings(geneId, closestBreakend.position());

                        exonsSkipped = abs(rnaExonData[EXON_RANK_MIN] - svPosExonData[EXON_RANK_MIN]);
                    }

                    rnaFusion.setTranscriptData(fs, closestTrans, closestBreakend, isViable, correctLocation, exonsSkipped);

                    LNX_LOGGER.debug("rnaFusion({}) {} closest breakend({}) distance({})",
                            rnaFusion.name(), isUpstream ? "up" :"down", closestBreakend.toString(), closestDistance);
                }
            }
        }

        setDnaFusionMatch(rnaFusion);
        rnaFusion.setFusionClusterChainInfo();
    }

    private void setDnaFusionMatch(final RnaFusionData rnaFusionData)
    {
        DnaRnaMatchType matchType = DnaRnaMatchType.NONE;

        final Transcript transUp = rnaFusionData.getMatchedfTranscripts()[FS_UPSTREAM];
        final Transcript transDown = rnaFusionData.getMatchedfTranscripts()[FS_DOWNSTREAM];

        for(final GeneFusion dnaFusion : mDnaFusions)
        {
            if(!dnaFusion.upstreamTrans().geneName().equals(rnaFusionData.GeneNames[FS_UPSTREAM])
            || !dnaFusion.downstreamTrans().geneName().equals(rnaFusionData.GeneNames[FS_DOWNSTREAM]))
            {
                continue;
            }

            matchType = DnaRnaMatchType.GENES;

            if(transUp != null && dnaFusion.upstreamTrans().gene().id() == transUp.gene().id()
            && transDown != null && dnaFusion.downstreamTrans().gene().id() == transDown.gene().id())
            {
                matchType = DnaRnaMatchType.SVS;
                break;
            }
        }

        if(matchType != DnaRnaMatchType.NONE)
        {
            rnaFusionData.setDnaFusionMatch(matchType, "");
            return;
        }

        for(Map.Entry<GeneFusion,String> entry : mDnaInvalidFusions.entrySet())
        {
            if(entry.getKey().name().equals(rnaFusionData.name()))
            {
                rnaFusionData.setDnaFusionMatch(DnaRnaMatchType.INVALID, entry.getValue());
            }
        }
    }

    private static int MAX_PROMOTOR_DISTANCE_UP = 100000;

    public boolean isTranscriptBreakendViableForRnaBoundary(final Transcript trans, boolean isUpstream, int breakendPosition,
            int rnaPosition, boolean exactRnaPosition)
    {
        // breakend must fall at or before the RNA boundary but not further upstream than the previous splice acceptor

        // if the RNA boundary is at or before the 2nd exon (which has the first splice acceptor), then the breakend can
        // be upstream as far the previous gene or 100K
        final TranscriptData transData = mGeneTransCache.getTranscriptData(trans.gene().StableId, trans.StableId);

        if (transData == null || transData.exons().isEmpty())
            return false;

        int strand = trans.gene().Strand;

        // first find the matching exon boundary for this RNA fusion boundary
        for (int i = 0; i < transData.exons().size(); ++i)
        {
            final ExonData exonData = transData.exons().get(i);
            final ExonData prevExonData = i > 0 ? transData.exons().get(i - 1) : null;
            final ExonData nextExonData = i < transData.exons().size() - 1 ? transData.exons().get(i + 1) : null;

            if (isUpstream)
            {
                // first check if at an exon boundary or before the start of the next exon and after the start of this one
                if(strand == 1)
                {
                    if ((rnaPosition == exonData.ExonEnd)
                    || (!exactRnaPosition && nextExonData != null && rnaPosition > exonData.ExonStart && rnaPosition < nextExonData.ExonStart))
                    {
                        // in which case check whether the breakend is before the next exon's splice acceptor
                        if (nextExonData != null)
                        {
                            return breakendPosition < nextExonData.ExonStart;
                        }

                        // can't take the last exon
                        return false;
                    }
                }
                else
                {
                    if ((rnaPosition == exonData.ExonStart)
                    || (!exactRnaPosition && prevExonData != null && rnaPosition < exonData.ExonEnd && rnaPosition > prevExonData.ExonEnd))
                    {
                        if(prevExonData != null)
                        {
                            return breakendPosition > prevExonData.ExonEnd;
                        }

                        return false;
                    }
                }
            }
            else
            {
                if((strand == 1 && rnaPosition <= exonData.ExonStart && exonData.ExonRank <= 2)
                || (strand == -1 && rnaPosition >= exonData.ExonEnd && exonData.ExonRank <= 2))
                {
                    int breakendDistance = abs(breakendPosition - rnaPosition);

                    if(breakendDistance > MAX_PROMOTOR_DISTANCE_UP || trans.hasNegativePrevSpliceAcceptorDistance())
                        return false;
                    else
                        return true;
                }

                if(strand == 1)
                {
                    if ((rnaPosition == exonData.ExonStart)
                    || (!exactRnaPosition && prevExonData != null && rnaPosition > prevExonData.ExonStart && rnaPosition < exonData.ExonStart))
                    {
                        if(prevExonData != null)
                        {
                            // after the previous exon's splice acceptor
                            return breakendPosition > prevExonData.ExonStart;
                        }

                        return false;
                    }
                }
                else
                {
                    if ((rnaPosition == exonData.ExonEnd)
                    || (!exactRnaPosition && nextExonData != null && rnaPosition < nextExonData.ExonEnd && rnaPosition > exonData.ExonEnd))
                    {
                        if(nextExonData != null)
                        {
                            // after the previous exon's splice acceptor
                            return breakendPosition < nextExonData.ExonStart;
                        }

                        return false;
                    }
                }
            }
        }

        return false;
    }

    private boolean isCandidateBetter(
            final GeneFusion currentFusion, final SvBreakend beCurrentStart, final SvBreakend beCurrentEnd,
            final GeneFusion candidateFusion, final SvBreakend beCandidateStart, final SvBreakend beCandidateEnd,
            final RnaFusionData rnaFusion, boolean currentFusionViable, boolean candidateFusionViable)
    {
        // if all else is equal, take a viable fusion over one that isn't
        if(beCurrentStart == beCandidateStart && beCurrentEnd == beCandidateEnd)
        {
            if (currentFusionViable != candidateFusionViable)
            {
                return candidateFusionViable;
            }

            if (currentFusion.phaseMatched() != candidateFusion.phaseMatched())
            {
                return candidateFusion.phaseMatched();
            }

            return false;
        }

        SvVarData currentStartSV = beCurrentStart.getSV();
        SvVarData currentEndSV = beCurrentEnd.getSV();
        SvVarData candidateStartSV = beCandidateStart.getSV();
        SvVarData candidateEndSV = beCandidateEnd.getSV();

        // give priority to same SV
        boolean currentSameSV = currentStartSV == currentEndSV;
        boolean candidateSameSV = candidateStartSV == candidateEndSV;

        if(currentSameSV != candidateSameSV)
            return candidateSameSV;

        // then whether chained
        boolean currentSameCluster = currentStartSV.getCluster() == currentEndSV.getCluster();
        boolean candidateSameCluster = candidateStartSV.getCluster() == candidateEndSV.getCluster();

        if(currentSameCluster != candidateSameCluster)
        {
            LNX_LOGGER.trace("current pair({} & {}) clusters({} & {}), candidate pair({} & {}) clusters({} & {})",
                    currentStartSV.id(), currentEndSV.id(), currentStartSV.getCluster().id(), currentEndSV.getCluster().id(),
                    candidateStartSV.id(), candidateEndSV.id(), candidateStartSV.getCluster().id(), candidateEndSV.getCluster().id());

            return candidateSameCluster;
        }

        if(currentSameCluster && candidateSameCluster)
        {
            // check whether one pair is in the same chain and the other not
            SvChain currentMatchingChain = currentStartSV.getCluster().findSameChainForSVs(currentStartSV, currentEndSV);

            SvChain candidateMatchingChain = candidateStartSV.getCluster().findSameChainForSVs(candidateStartSV, candidateEndSV);

            LNX_LOGGER.trace("current pair({} & {}) clusters({} chain={}), candidate pair({} & {}) clusters({} chain={})",
                    currentStartSV.id(), currentEndSV.id(), currentStartSV.getCluster().id(),
                    currentMatchingChain != null ? currentMatchingChain.id() : "diff",
                    candidateStartSV.id(), candidateEndSV.id(), candidateStartSV.getCluster().id(),
                    candidateMatchingChain != null ? candidateMatchingChain.id() : "diff");

            if(currentMatchingChain != null && candidateMatchingChain == null)
                return false;
            if(currentMatchingChain == null && candidateMatchingChain != null)
                return true;
        }

        // otherwise revert to whichever positions are closest to the RNA breakends

        // lastly the nearest to the RNA positions
        double currentPosDiff = (abs(rnaFusion.Positions[FS_UPSTREAM] - beCurrentStart.position())
                + abs(rnaFusion.Positions[FS_DOWNSTREAM] - beCurrentEnd.position())) * 0.5;

        double candidatePosDiff = (abs(rnaFusion.Positions[FS_UPSTREAM] - beCandidateStart.position())
                + abs(rnaFusion.Positions[FS_DOWNSTREAM] - beCandidateEnd.position())) * 0.5;

        return candidatePosDiff < currentPosDiff;
    }

    private boolean isViableBreakend(final SvBreakend breakend, int rnaPosition, byte geneStrand, boolean isUpstream)
    {
        boolean requireHigherBreakendPos = isUpstream ? (geneStrand == 1) : (geneStrand == -1);

        int position = breakend.position();

        int offsetMargin = breakend.usesStart() ?
                breakend.getSV().getSvData().startHomologySequence().length() : breakend.getSV().getSvData().endHomologySequence().length();

        offsetMargin = offsetMargin / 2 + 1;

        // the interval offset could be used in place of half the homology but interpretation of the GRIDSS value needs to be understood first
        /*
        int offsetMargin = requireHigherBreakendPos ?
                (breakend.usesStart() ? breakend.getSV().getSvData().startIntervalOffsetEnd() : breakend.getSV().getSvData().endIntervalOffsetEnd())
                : (breakend.usesStart() ? breakend.getSV().getSvData().startIntervalOffsetStart() : breakend.getSV().getSvData().endIntervalOffsetStart());
        */

        if(requireHigherBreakendPos)
        {
            if(breakend.orientation() != 1)
                return false;

            // factor in any uncertainty around the precise breakend, eg from homology
            return (position + offsetMargin >= rnaPosition);
        }
        else
        {
            if(breakend.orientation() != -1)
                return false;

            return (position - offsetMargin <= rnaPosition);
        }
    }

    private void setRnaFusionData(final RnaFusionData rnaFusion)
    {
        // find transcripts which match the RNA positions
        Map<String,int[]> transPhasesUp = Maps.newHashMap();
        Map<String,int[]> transPhasesDown = Maps.newHashMap();

        boolean isExactRnaExon = rnaFusion.matchesKnownSpliceSite();

        for(int fs = FS_UPSTREAM; fs <= 1; ++fs)
        {
            boolean isUpstream = (fs == 0);
            final String geneName = rnaFusion.GeneNames[fs];
            int rnaPosition = rnaFusion.Positions[fs];

            Map<String,int[]> transPhases = isUpstream ? transPhasesUp : transPhasesDown;

            EnsemblGeneData geneData = null;

            if(rnaFusion.GeneIds[fs].isEmpty())
            {
                geneData = mGeneTransCache.getGeneDataByName(rnaFusion.GeneNames[fs]);

                if(geneData == null)
                {
                    LNX_LOGGER.warn("sample({}) rnaFusion({}) {} gene not found", mSampleId, rnaFusion.name(), isUpstream ? "up" : "down");
                    rnaFusion.setValid(false);
                    return;
                }

                rnaFusion.GeneIds[fs] = geneData.GeneId;
            }
            else
            {
                geneData = mGeneTransCache.getGeneDataById(rnaFusion.GeneIds[fs]);
            }

            rnaFusion.Strands[fs] = geneData.Strand;

            if (rnaFusion.JunctionTypes[fs] == RnaJunctionType.KNOWN)
            {
                // check that the RNA position is within the bounds of the gene before proceeding
                int upPosLimit = geneData.GeneEnd;
                int downPosLimit = geneData.GeneStart;

                if(!isUpstream)
                {
                    if(geneData.Strand == 1)
                        downPosLimit -= PRE_GENE_PROMOTOR_DISTANCE;
                    else
                        upPosLimit += PRE_GENE_PROMOTOR_DISTANCE;
                }

                if(rnaPosition < downPosLimit || rnaPosition > upPosLimit)
                {
                    LNX_LOGGER.warn("sample({}) rnaFusion({}) {} position({}) outside geneBounds({} -> {})",
                            mSampleId, rnaFusion.name(), isUpstream ? "upstream" : "downstream", rnaPosition, downPosLimit, upPosLimit);
                    rnaFusion.setValid(false);
                    return;
                }
            }

            List<TranscriptData> transDataList = mGeneTransCache.getTranscripts(geneData.GeneId);

            if(transDataList != null)
            {
                for(final TranscriptData transData : transDataList)
                {
                    // skip any downstream gene 3'UTR or non-coding transcripts for this position
                    if(!isUpstream)
                    {
                        if (transData.CodingStart == null)
                            continue;

                        if (transData.Strand == 1 && transData.CodingEnd != null && rnaPosition >= transData.CodingEnd)
                            continue;
                        else if (transData.Strand == -1 && transData.CodingStart != null && rnaPosition <= transData.CodingStart)
                            continue;
                    }

                    int[] exonMatchData = findExonMatch(transData.exons(), transData.Strand, rnaPosition);

                    if(exonMatchData[EXON_EXACT_MATCH] > 0)
                        rnaFusion.getExactMatchTransIds()[fs].add(transData.TransName);

                    if(exonMatchData[EXON_FOUND] > 0)
                    {
                        transPhases.put(transData.TransName, exonMatchData);
                    }
                }
            }

            LNX_LOGGER.debug("rnaFusion({}) juncType({}) {} position({}) matched {} transcripts",
                    rnaFusion.name(), rnaFusion.JunctionTypes[fs], streamStr(fs), rnaPosition, transPhases.size());
        }

        if(!transPhasesUp.isEmpty() && !transPhasesDown.isEmpty())
        {
            for (Map.Entry<String, int[]> entryUp : transPhasesUp.entrySet())
            {
                final String transIdUp = entryUp.getKey();
                final int[] exonDataUp = entryUp.getValue();

                if(isExactRnaExon && !rnaFusion.getExactMatchTransIds()[FS_UPSTREAM].contains(transIdUp))
                    continue;

                for (Map.Entry<String, int[]> entryDown : transPhasesDown.entrySet())
                {
                    final String transIdDown = entryDown.getKey();
                    final int[] exonDataDown = entryDown.getValue();

                    if(isExactRnaExon && !rnaFusion.getExactMatchTransIds()[FS_DOWNSTREAM].contains(transIdDown))
                        continue;

                    boolean phaseMatched = exonDataUp[EXON_PHASE] == exonDataDown[EXON_PHASE];

                    if (phaseMatched && !rnaFusion.hasRnaPhasedFusion())
                    {
                        LNX_LOGGER.debug("rnaFusion({}) juncTypes(up={} down={}) transUp({}) transDown({} phase({}) matched({})",
                                rnaFusion.name(), rnaFusion.JunctionTypes[FS_UPSTREAM], rnaFusion.JunctionTypes[FS_DOWNSTREAM],
                                transIdUp, transIdDown, exonDataUp[EXON_PHASE], phaseMatched);

                        rnaFusion.setRnaPhasedFusionData(transIdUp, transIdDown);
                        rnaFusion.setExonData(FS_UPSTREAM, exonDataUp[EXON_RANK], exonDataUp[EXON_PHASE]);
                        rnaFusion.setExonData(FS_DOWNSTREAM, exonDataDown[EXON_RANK], exonDataDown[EXON_PHASE]);
                    }
                }
            }
        }

        KnownFusionData refFusionData = mFusionFinder.getKnownFusionData();

        if(refFusionData == null)
        {
            rnaFusion.setKnownType(REPORTABLE_TYPE_NONE);
            return;
        }

        for(final String[] genePair : refFusionData.knownPairs())
        {
            if (genePair[FIVE_GENE].equals(rnaFusion.GeneNames[FS_UPSTREAM]) && genePair[THREE_GENE].equals(rnaFusion.GeneNames[FS_DOWNSTREAM]))
            {
                rnaFusion.setKnownType(REPORTABLE_TYPE_KNOWN);
                return;
            }
        }

        boolean fivePrimeProm = refFusionData.hasPromiscuousFiveGene(rnaFusion.GeneNames[FS_UPSTREAM]);
        boolean threePrimeProm = refFusionData.hasPromiscuousThreeGene(rnaFusion.GeneNames[FS_DOWNSTREAM]);

        if(fivePrimeProm && threePrimeProm)
        {
            rnaFusion.setKnownType(REPORTABLE_TYPE_BOTH_PROM);
        }
        else if(fivePrimeProm)
        {
            rnaFusion.setKnownType(REPORTABLE_TYPE_5P_PROM);
        }
        else if(threePrimeProm)
        {
            rnaFusion.setKnownType(REPORTABLE_TYPE_3P_PROM);
        }
        else
        {
            rnaFusion.setKnownType(REPORTABLE_TYPE_NONE);
        }
    }

    private static final int EXON_FOUND = 0;
    private static final int EXON_RANK = 1;
    private static final int EXON_PHASE = 2;
    private static final int EXON_EXACT_MATCH = 3;

    private int[] findExonMatch(final List<ExonData> exonDataList, int strand, int position)
    {
        int[] exonMatch = new int[EXON_EXACT_MATCH+1];

        for (int i = 0; i < exonDataList.size(); ++i)
        {
            final ExonData transExonData = exonDataList.get(i);
            final ExonData nextTransExonData = i < exonDataList.size() - 1 ? exonDataList.get(i + 1) : null;

            if (position == transExonData.ExonEnd || position == transExonData.ExonStart)
            {
                // skip matches on the last exon
                if(i == 0 && strand == -1 && position == transExonData.ExonStart)
                    return exonMatch;
                else if(i == exonDataList.size() - 1 && strand == 1 && position == transExonData.ExonEnd)
                    return exonMatch;

                // position exactly matches the bounds of an exon
                exonMatch[EXON_FOUND] = 1;
                exonMatch[EXON_EXACT_MATCH] = 1;
                exonMatch[EXON_RANK] = transExonData.ExonRank;

                if ((strand == 1) == (position == transExonData.ExonStart))
                {
                    exonMatch[EXON_PHASE] = transExonData.ExonPhase;
                }
                else
                {
                    exonMatch[EXON_PHASE] = transExonData.ExonPhaseEnd;
                }
                break;
            }

            if (position > transExonData.ExonStart && position < transExonData.ExonEnd)
            {
                // position is within the bounds of an exon
                exonMatch[EXON_FOUND] = 1;
                exonMatch[EXON_RANK] = transExonData.ExonRank;
                exonMatch[EXON_PHASE] = transExonData.ExonPhase;
                break;
            }

            if (nextTransExonData != null && position > transExonData.ExonEnd && position < nextTransExonData.ExonStart)
            {
                exonMatch[EXON_FOUND] = 1;
                if (strand == 1)
                {
                    exonMatch[EXON_RANK] = transExonData.ExonRank;
                    exonMatch[EXON_PHASE] = transExonData.ExonPhaseEnd;
                }
                else
                {
                    exonMatch[EXON_RANK] = nextTransExonData.ExonRank;
                    exonMatch[EXON_PHASE] = nextTransExonData.ExonPhaseEnd;
                }

                break;
            }
        }

        return exonMatch;
    }

    public boolean loadSampleRnaData(final String source, final String filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));

            final Map<String,Integer> fieldIndexMap = createFieldsIndexMap(lines.get(0), getRnaSourceDelimiter(source));

            lines.remove(0);

            String currentSampleId = "";
            List<RnaFusionData> rnaDataList = Lists.newArrayList();
            int recordCount = 0;

            for(String data : lines)
            {
                RnaFusionData rnaData = RnaFusionData.from(source, recordCount, data, fieldIndexMap);
                ++recordCount;

                if(currentSampleId.isEmpty() || !currentSampleId.equals(rnaData.SampleId))
                {
                    currentSampleId = rnaData.SampleId;
                    rnaDataList = Lists.newArrayList();
                    mSampleRnaData.put(currentSampleId, rnaDataList);
                }

                rnaDataList.add(rnaData);
            }
        }
        catch(IOException e)
        {
            LNX_LOGGER.warn("failed to load sample RNA data file({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

    public void close() { mWriter.close(); }

}
