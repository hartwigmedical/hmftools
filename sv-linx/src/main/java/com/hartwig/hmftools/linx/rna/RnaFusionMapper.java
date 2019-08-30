package com.hartwig.hmftools.linx.rna;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_3P_PROM;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_5P_PROM;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_BOTH_PROM;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_KNOWN;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_NONE;
import static com.hartwig.hmftools.linx.fusion.KnownFusionData.FIVE_GENE;
import static com.hartwig.hmftools.linx.fusion.KnownFusionData.THREE_GENE;
import static com.hartwig.hmftools.linx.rna.RnaFusionData.RNA_SPLICE_TYPE_ONLY_REF;
import static com.hartwig.hmftools.linx.rna.RnaFusionData.RNA_SPLICE_TYPE_UNKONWN;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.isSpecificSV;
import static com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection.EXON_RANK_MIN;
import static com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection.PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.linx.fusion.FusionFinder.checkFusionLogic;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.linx.fusion.FusionParameters;
import com.hartwig.hmftools.linx.fusion.KnownFusionData;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;
import com.hartwig.hmftools.linx.fusion.FusionFinder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class RnaFusionMapper
{
    private String mSampleId;
    private String mOutputDir;

    private FusionFinder mFusionFinder;
    private FusionParameters mFusionParams;
    private SvGeneTranscriptCollection mGeneTransCollection;
    private Map<String, List<RnaFusionData>> mSampleRnaData;
    private final List<GeneFusion> mDnaFusions;

    private BufferedWriter mWriter;

    private static final Logger LOGGER = LogManager.getLogger(RnaFusionMapper.class);

    public RnaFusionMapper(SvGeneTranscriptCollection geneTransCollection, FusionFinder fusionFinder, final List<GeneFusion> dnaFusions)
    {
        mSampleRnaData = Maps.newHashMap();
        mWriter = null;
        mFusionFinder = fusionFinder;
        mGeneTransCollection = geneTransCollection;
        mDnaFusions = dnaFusions;

        mFusionParams = new FusionParameters();
        mFusionParams.RequirePhaseMatch = false;
        mFusionParams.AllowExonSkipping = false;
    }

    public void setOutputDir(final String outputDir)
    {
        mOutputDir = outputDir;
    }

    public final Map<String, List<RnaFusionData>> getSampleRnaData() { return mSampleRnaData; }
    public final List<RnaFusionData> getSampleRnaData(final String sampleId) { return mSampleRnaData.get(sampleId); }

    public void assessRnaFusions(final String sampleId, Map<String, List<SvBreakend>> chrBreakendMap)
    {
        mSampleId = sampleId;

        final List<RnaFusionData> rnaFusionList = mSampleRnaData.get(mSampleId);

        if (rnaFusionList == null || rnaFusionList.isEmpty())
            return;

        LOGGER.debug("assessing {} RNA fusions", rnaFusionList.size());

        for (final RnaFusionData rnaFusion : rnaFusionList)
        {
            setRnaFusionData(rnaFusion);

            if(!rnaFusion.isValid())
                continue;

            annotateRnaFusions(rnaFusion, chrBreakendMap);

            writeRnaMatchData(mSampleId, rnaFusion);
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

        boolean isExactRnaExon = rnaFusion.SpliceType.equals(RNA_SPLICE_TYPE_ONLY_REF);

        for(int i = 0; i <= 1 ; ++i)
        {
            boolean isUpstream = (i == 0);
            String chromosome = isUpstream ? rnaFusion.ChrUp : rnaFusion.ChrDown;
            long rnaPosition = isUpstream ? rnaFusion.PositionUp : rnaFusion.PositionDown;
            byte geneStrand = isUpstream ? rnaFusion.StrandUp : rnaFusion.StrandDown;
            List<SvBreakend> viableBreakends = isUpstream ? viableUpBreakends : viableDownBreakends;
            List<SvBreakend> nearBreakends = isUpstream ? nearUpBreakends : nearDownBreakends;
            List<SvBreakend> genicBreakends = isUpstream ? genicUpBreakends : genicDownBreakends;
            List<Transcript> viableTranscripts = isUpstream ? viableUpTranscripts : viableDownTranscripts;
            List<Transcript> nearTranscripts = isUpstream ? nearUpTranscripts : nearDownTranscripts;
            List<Transcript> genicTranscripts = isUpstream ? genicUpTranscripts : genicDownTranscripts;
            String geneName = isUpstream ? rnaFusion.GeneUp : rnaFusion.GeneDown;

            final List<SvBreakend> breakendList = chrBreakendMap.get(chromosome);

            if(breakendList == null)
                continue;

            for(final SvBreakend breakend : breakendList)
            {
                final SvVarData var = breakend.getSV();

                if(var.isInferredSgl())
                    continue;

                isSpecificSV(var);

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
                        if(!rnaFusion.getExactMatchTransIds(isUpstream).contains(trans.TransId))
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

        LOGGER.debug("rna fusion({}) breakend matches: upstream(viable={} near={} genic={}) downstream(viable={} near={} genic={})",
                rnaFusion.Name, viableUpBreakends.size(), nearUpBreakends.size(), genicUpBreakends.size(),
                viableDownBreakends.size(), nearDownBreakends.size(), genicDownBreakends.size());

        // run them through fusion logic (ie a pair of breakend lists), but don't require phase matching
        if(!viableUpBreakends.isEmpty() && !viableDownBreakends.isEmpty())
        {
            GeneFusion topCandidateFusion = null;
            SvBreakend topUpBreakend = null;
            SvBreakend topDownBreakend = null;

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

                    // form one any way but mark it as not meeting standard fusion rules
                    if(possibleFusion == null)
                    {
                        possibleFusion = new GeneFusion(upTrans, downTrans, false, false);
                    }
                    else if(!matchesDnaFusion(possibleFusion))
                    {
                        possibleFusion.setViable(false);
                    }

                    if (topCandidateFusion == null
                    || isCandidateBetter(topCandidateFusion, topUpBreakend, topDownBreakend, possibleFusion, upBreakend, downBreakend, rnaFusion))
                    {
                        topCandidateFusion = possibleFusion;
                        topUpBreakend = upBreakend;
                        topDownBreakend = downBreakend;

                        LOGGER.debug("rnaFusion({}) first pair({} & {})", rnaFusion.Name, upBreakend.toString(), downBreakend.toString());
                    }
                }
            }

            if(topCandidateFusion != null)
            {
                rnaFusion.setTranscriptData(
                        true, topCandidateFusion.upstreamTrans(), topUpBreakend,
                        true, true,  0);

                rnaFusion.setTranscriptData(
                        false, topCandidateFusion.downstreamTrans(), topDownBreakend,
                        true, true,0);

                rnaFusion.setViableFusion(topCandidateFusion.viable(), topCandidateFusion.phaseMatched());
            }
        }
        else
        {
            // select the closest breakend's transcript
            for(int i = 0; i <= 1 ; ++i)
            {
                boolean isUpstream = (i == 0);
                long rnaPosition = isUpstream ? rnaFusion.PositionUp : rnaFusion.PositionDown;

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
                long closestDistance = 0;

                for (int j = 0; j < transcriptList.size(); ++j)
                {
                    final Transcript trans = transcriptList.get(j);
                    final SvBreakend breakend = breakendList.get(j);

                    long distance = abs(rnaPosition - trans.svPosition());
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
                        final int rnaExonData[] = mGeneTransCollection.getExonRankings(geneId, rnaPosition);
                        final int svPosExonData[] = mGeneTransCollection.getExonRankings(geneId, closestBreakend.position());

                        exonsSkipped = abs(rnaExonData[EXON_RANK_MIN] - svPosExonData[EXON_RANK_MIN]);
                    }

                    rnaFusion.setTranscriptData(isUpstream, closestTrans, closestBreakend, isViable, correctLocation, exonsSkipped);

                    LOGGER.debug("rnaFusion({}) {} closest breakend({}) distance({})",
                            rnaFusion.Name, isUpstream ? "up" :"down", closestBreakend.toString(), closestDistance);
                }
            }
        }

        rnaFusion.setFusionClusterChainInfo();
    }

    private boolean matchesDnaFusion(final GeneFusion fusion)
    {
        for(final GeneFusion dnaFusion : mDnaFusions)
        {
            if(dnaFusion.upstreamTrans().gene().id() == fusion.upstreamTrans().gene().id()
            && dnaFusion.downstreamTrans().gene().id() == fusion.downstreamTrans().gene().id())
            {
                return true;
            }
        }

        return false;
    }

    private static int MAX_PROMOTOR_DISTANCE_UP = 100000;

    public boolean isTranscriptBreakendViableForRnaBoundary(final Transcript trans, boolean isUpstream, long breakendPosition,
            long rnaPosition, boolean exactRnaPosition)
    {
        // breakend must fall at or before the RNA boundary but not further upstream than the previous splice acceptor

        // if the RNA boundary is at or before the 2nd exon (which has the first splice acceptor), then the breakend can
        // be upstream as far the previous gene or 100K
        final TranscriptData transData = mGeneTransCollection.getTranscriptData(trans.gene().StableId, trans.StableId);

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
                    int distanceUp = trans.prevSpliceAcceptorDistance();
                    long breakendDistance = abs(breakendPosition - rnaPosition);

                    if(breakendDistance > MAX_PROMOTOR_DISTANCE_UP || distanceUp < 0)
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

    private boolean isCandidateBetter(final GeneFusion currentFusion, final SvBreakend beCurrentStart, final SvBreakend beCurrentEnd,
            final GeneFusion candidateFusion, final SvBreakend beCandidateStart, final SvBreakend beCandidateEnd, final RnaFusionData rnaFusion)
    {
        // if all else is equal, take a viable fusion over one that isn't
        if(beCurrentStart == beCandidateStart && beCurrentEnd == beCandidateEnd)
        {
            if (currentFusion.viable() != candidateFusion.viable())
            {
                return candidateFusion.viable();
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
            LOGGER.debug("current pair({} & {}) clusters({} & {}), candidate pair({} & {}) clusters({} & {})",
                    currentStartSV.id(), currentEndSV.id(), currentStartSV.getCluster().id(), currentEndSV.getCluster().id(),
                    candidateStartSV.id(), candidateEndSV.id(), candidateStartSV.getCluster().id(), candidateEndSV.getCluster().id());

            return candidateSameCluster;
        }

        if(currentSameCluster && candidateSameCluster)
        {
            // check whether one pair is in the same chain and the other not
            SvChain currentMatchingChain = currentStartSV.getCluster().findSameChainForSVs(currentStartSV, currentEndSV);

            SvChain candidateMatchingChain = candidateStartSV.getCluster().findSameChainForSVs(candidateStartSV, candidateEndSV);

            LOGGER.debug("current pair({} & {}) clusters({} chain={}), candidate pair({} & {}) clusters({} chain={})",
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
        double currentPosDiff = (abs(rnaFusion.PositionUp - beCurrentStart.position()) + abs(rnaFusion.PositionDown - beCurrentEnd.position())) * 0.5;
        double candidatePosDiff = (abs(rnaFusion.PositionUp - beCandidateStart.position()) + abs(rnaFusion.PositionDown - beCandidateEnd.position())) * 0.5;

        return candidatePosDiff < currentPosDiff;
    }

    private boolean isViableBreakend(final SvBreakend breakend, long rnaPosition, byte geneStrand, boolean isUpstream)
    {
        boolean requireHigherBreakendPos = isUpstream ? (geneStrand == 1) : (geneStrand == -1);

        long position = breakend.position();

        if(requireHigherBreakendPos)
        {
            // factor in any uncertainty around the precise breakend, eg from homology
            position += breakend.usesStart() ? breakend.getSV().getSvData().startIntervalOffsetEnd() : breakend.getSV().getSvData().endIntervalOffsetEnd();

            return (position >= rnaPosition);
        }
        else
        {
            position += breakend.usesStart() ? breakend.getSV().getSvData().startIntervalOffsetStart() : breakend.getSV().getSvData().endIntervalOffsetStart();

            return (position <= rnaPosition);
        }
    }

    private void setRnaFusionData(final RnaFusionData rnaFusion)
    {
        // find transcripts which match the RNA positions
        Map<String,Integer> transPhasesUp = Maps.newHashMap();
        Map<String,Integer> transPhasesDown = Maps.newHashMap();

        for(int i = 0; i <= 1; ++i)
        {
            boolean isUpstream = (i == 0);
            final String geneName = isUpstream ? rnaFusion.GeneUp : rnaFusion.GeneDown;
            long rnaPosition = isUpstream ? rnaFusion.PositionUp : rnaFusion.PositionDown;
            Map<String,Integer> transPhases = isUpstream ? transPhasesUp : transPhasesDown;

            EnsemblGeneData geneData = mGeneTransCollection.getGeneDataByName(geneName);

            if (geneData == null)
            {
                LOGGER.warn("sample({}) rnaFusion({}) {} gene not found", mSampleId, rnaFusion.Name, isUpstream ? "up" : "down");
                rnaFusion.setValid(false);
                return;
            }

            rnaFusion.setGeneId(geneData.GeneId, isUpstream);

            if (rnaFusion.SpliceType.equals(RNA_SPLICE_TYPE_UNKONWN))
            {
                // check that the RNA position is within the bounds of the gene before proceeding
                long upPosLimit = geneData.GeneEnd;
                long downPosLimit = geneData.GeneStart;

                if(!isUpstream)
                {
                    if(geneData.Strand == 1)
                        downPosLimit -= PRE_GENE_PROMOTOR_DISTANCE;
                    else
                        upPosLimit += PRE_GENE_PROMOTOR_DISTANCE;
                }

                if(rnaPosition < downPosLimit || rnaPosition > upPosLimit)
                {
                    LOGGER.warn("sample({}) rnaFusion({}) {} position({}) outside geneBounds({} -> {})",
                            mSampleId, rnaFusion.Name, isUpstream ? "upstream" : "downstream", rnaPosition, downPosLimit, upPosLimit);
                    rnaFusion.setValid(false);
                    return;
                }
            }

            List<TranscriptData> transDataList = mGeneTransCollection.getTranscripts(geneData.GeneId);

            if(transDataList != null)
            {
                for(final TranscriptData transData : transDataList)
                {
                    int[] exonMatchData = findExonMatch(transData.exons(), transData.Strand, rnaPosition);

                    if(exonMatchData[EXON_EXACT_MATCH] > 0)
                        rnaFusion.getExactMatchTransIds(isUpstream).add(transData.TransId);

                    if(exonMatchData[EXON_FOUND] > 0)
                    {
                        transPhases.put(transData.TransName, exonMatchData[EXON_PHASE]);

                        if(transData.IsCanonical || rnaFusion.exonRank(isUpstream) == 0)
                        {
                            rnaFusion.setExonData(isUpstream, exonMatchData[EXON_RANK],exonMatchData[EXON_PHASE]);
                        }
                    }
                }
            }

            LOGGER.debug("rnaFusion({}) type({}) {} position({}) matched {} transcripts",
                    rnaFusion.Name, rnaFusion.SpliceType, isUpstream ? "upstream" : "downstream", rnaPosition, transPhases.size());

            // int[] transUpExonData = mGeneTransCollection.getExonRankings(geneData.GeneId, rnaPosition);
            // rnaFusion.setExonRank(isUpstream, transUpExonData[EXON_RANK_MIN], transUpExonData[EXON_RANK_MAX]);
        }

        if(!transPhasesUp.isEmpty() && !transPhasesDown.isEmpty())
        {
            for (Map.Entry<String, Integer> entryUp : transPhasesUp.entrySet())
            {
                final String transIdUp = entryUp.getKey();
                int exonPhaseUp = entryUp.getValue();

                for (Map.Entry<String, Integer> entryDown : transPhasesDown.entrySet())
                {
                    final String transIdDown = entryDown.getKey();
                    int exonPhaseDown = entryDown.getValue();

                    boolean phaseMatched = exonPhaseUp == exonPhaseDown;

                    if (phaseMatched && !rnaFusion.hasRnaPhasedFusion())
                    {
                        LOGGER.debug("rnaFusion({}) type({}) transUp({}) transDown({} phase({}) matched({})",
                                rnaFusion.Name, rnaFusion.SpliceType, transIdUp, transIdDown, exonPhaseUp, phaseMatched);

                        rnaFusion.setRnaPhasedFusionData(transIdUp, transIdDown);
                    }
                }
            }
        }

        KnownFusionData refFusionData = mFusionFinder.getKnownFusionDatal();

        if(refFusionData == null)
        {
            rnaFusion.setKnownType(REPORTABLE_TYPE_NONE);
            return;
        }

        for(final String[] genePair : refFusionData.knownPairs())
        {
            if (genePair[FIVE_GENE].equals(rnaFusion.GeneUp) && genePair[THREE_GENE].equals(rnaFusion.GeneDown))
            {
                rnaFusion.setKnownType(REPORTABLE_TYPE_KNOWN);
                return;
            }
        }

        boolean fivePrimeProm = refFusionData.hasPromiscuousFiveGene(rnaFusion.GeneUp);
        boolean threePrimeProm = refFusionData.hasPromiscuousThreeGene(rnaFusion.GeneDown);

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

    private int[] findExonMatch(final List<ExonData> exonDataList, int strand, long position)
    {
        int[] exonMatch = new int[EXON_EXACT_MATCH+1];

        for (int i = 0; i < exonDataList.size(); ++i)
        {
            final ExonData transExonData = exonDataList.get(i);
            final ExonData nextTransExonData = i < exonDataList.size() - 1 ? exonDataList.get(i + 1) : null;

            if (position == transExonData.ExonEnd || position == transExonData.ExonStart)
            {
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
                    exonMatch[EXON_PHASE] = transExonData.ExonPhase;
                }
                else
                {
                    exonMatch[EXON_RANK] = nextTransExonData.ExonRank;
                    exonMatch[EXON_PHASE] = nextTransExonData.ExonPhase;
                }

                break;
            }
        }

        return exonMatch;
    }

    public void writeRnaMatchData(final String sampleId, final RnaFusionData rnaFusion)
    {
        try
        {
            if(mWriter == null)
            {
                String outputFilename = mOutputDir;

                outputFilename += "SVA_RNA_DATA.csv";

                mWriter = createBufferedWriter(outputFilename, false);

                mWriter.write("SampleId,FusionName,GeneNameUp,GeneNameDown,ViableFusion,PhaseMatched,KnownType");

                for(int se = SE_START; se <= SE_END; ++se)
                {
                    String upDown = se == SE_START ? "Up" : "Down";

                    String fieldsStr = ",SvId" + upDown;
                    fieldsStr += ",Chr" + upDown;
                    fieldsStr += ",Pos" + upDown;
                    fieldsStr += ",RnaPos" + upDown;
                    fieldsStr += ",Orient" + upDown;
                    fieldsStr += ",Strand" + upDown;
                    fieldsStr += ",Type" + upDown;
                    fieldsStr += ",ClusterInfo" + upDown;
                    fieldsStr += ",TransViable" + upDown;
                    fieldsStr += ",TransValidLoc" + upDown;
                    fieldsStr += ",TransId" + upDown;
                    fieldsStr += ",ExonsSkipped" + upDown;
                    fieldsStr += ",RegionType" + upDown;
                    fieldsStr += ",CodingType" + upDown;
                    fieldsStr += ",Exon" + upDown;
                    fieldsStr += ",Disruptive" + upDown;
                    fieldsStr += ",DistancePrev" + upDown;
                    fieldsStr += ",RnaExonRank" + upDown;
                    fieldsStr += ",RnaExonPhase" + upDown;
                    fieldsStr += ",RnaExonMatch" + upDown;
                    mWriter.write(fieldsStr);
                }

                // mWriter.write(",SvIdUp,ChrUp,PosUp,RnaPosUp,OrientUp,StrandUp,TypeUp,ClusterInfoUp");
                // mWriter.write(",TransViableUp,TransValidLocUp,TransIdUp,ExonsSkippedUp,RegionTypeUp,CodingTypeUp,ExonUp,DisruptiveUp,DistancePrevUp");

                // mWriter.write(",SvIdDown,ChrDown,PosDown,RnaPosDown,OrientDown,StrandDown,TypeDown,ClusterInfoDown");
                // mWriter.write(",TransViableDown,TransValidLocDown,TransIdDown,ExonsSkippedDown,RegionTypeDown,CodingTypeDown,ExonDown,DisruptiveDown,DistancePrevDown");

                mWriter.write(",ChainInfo,JunctionReadCount,SpanningFragCount,SpliceType");
                mWriter.write(",RnaPhaseMatched,RnaTransIdUp,RnaTransIdDown");

                mWriter.newLine();
            }

            BufferedWriter writer = mWriter;

            writer.write(String.format("%s,%s,%s,%s,%s,%s,%s",
                    sampleId, rnaFusion.Name, rnaFusion.GeneUp, rnaFusion.GeneDown,
                    rnaFusion.isViableFusion(), rnaFusion.isPhaseMatchedFusion(), rnaFusion.getKnownType()));

            for(int se = SE_START; se <= SE_END; ++se)
            {
                boolean isUpstream = (se == SE_START);
                final Transcript trans = rnaFusion.getTrans(isUpstream);

                if(trans != null)
                {
                    final GeneAnnotation gene = trans.gene();

                    writer.write(String.format(",%d,%s,%d,%d,%d,%d,%s,%s",
                            gene.id(), gene.chromosome(), gene.position(),
                            isUpstream ? rnaFusion.PositionUp : rnaFusion.PositionDown,
                            gene.orientation(), gene.Strand, gene.type(),
                            rnaFusion.getClusterInfo(isUpstream)));

                    writer.write(String.format(",%s,%s,%s,%d,%s,%s,%d,%s,%d",
                            rnaFusion.isTransViable(isUpstream), rnaFusion.isTransCorrectLocation(isUpstream),
                            trans.StableId, rnaFusion.getExonsSkipped(isUpstream),
                            trans.regionType(), trans.codingType(),
                            isUpstream ? trans.ExonUpstream : trans.ExonDownstream, trans.isDisruptive(), trans.prevSpliceAcceptorDistance()));
                }
                else
                {
                    writer.write(String.format(",%s,%s,%d,%d,%d,%d,%s,%s",
                            "", isUpstream ? rnaFusion.ChrUp : rnaFusion.ChrDown, 0,
                            isUpstream ? rnaFusion.PositionUp : rnaFusion.PositionDown, 0,
                            isUpstream ? rnaFusion.StrandUp : rnaFusion.StrandDown, "", ""));

                    writer.write(String.format(",%s,%s,,,,,,,",
                            rnaFusion.isTransViable(isUpstream), rnaFusion.isTransCorrectLocation(isUpstream)));
                }

                boolean hasExonData = rnaFusion.exonRank(isUpstream) > 0;
                writer.write(String.format(",%s,%s,%d",
                        hasExonData ? rnaFusion.exonRank(isUpstream) :"",
                        hasExonData ? rnaFusion.exonPhase(isUpstream) : "",
                        rnaFusion.getExactMatchTransIds(isUpstream).size()));
            }

            writer.write(String.format(",%s,%d,%d,%s",
                    !rnaFusion.getChainInfo().isEmpty() ? rnaFusion.getChainInfo() : "0;0",
                    rnaFusion.JunctionReadCount, rnaFusion.SpanningFragCount, rnaFusion.SpliceType));

            writer.write(String.format(",%s,%s,%s", rnaFusion.hasRnaPhasedFusion(),
                    rnaFusion.getRnaPhasedFusionTransName(true), rnaFusion.getRnaPhasedFusionTransName(false)));

            writer.newLine();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing RNA match data: {}", e.toString());
        }
    }

    private static int COL_SAMPLEID = 0;
    private static int COL_NAME = 1;
    private static int COL_JUNCT_RC = 2;
    private static int COL_SPAN_RC = 3;
    private static int COL_SPLICE = 4;
    private static int COL_GENE_UP = 5;
    private static int COL_CHR_UP = 7;
    private static int COL_POS_UP = 8;
    private static int COL_STRAND_UP = 9;
    private static int COL_GENE_DOWN = 10;
    private static int COL_CHR_DOWN = 12;
    private static int COL_POS_DOWN = 13;
    private static int COL_STRAND_DOWN = 14;

    public boolean loadSampleRnaData(final String filename)
    {
        if (filename.isEmpty() || !Files.exists(Paths.get(filename)))
            return false;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("empty RNA data file({})", filename);
                return false;
            }

            line = fileReader.readLine(); // skip header

            String currentSampleId = "";
            List<RnaFusionData> rnaDataList = Lists.newArrayList();

            while (line != null)
            {
                // parse CSV data
                String[] items = line.split(",");

                // check if still on the same variant
                final String sampleId = items[COL_SAMPLEID];

                if(currentSampleId.isEmpty() || !currentSampleId.equals(sampleId))
                {
                    currentSampleId = sampleId;
                    rnaDataList = Lists.newArrayList();
                    mSampleRnaData.put(currentSampleId, rnaDataList);
                }

                // check that gene names match Ensembl
                String geneUp = items[COL_GENE_UP];
                String geneDown = items[COL_GENE_DOWN];

                geneUp = checkAlternateGeneName(geneUp);
                geneDown = checkAlternateGeneName(geneDown);

                RnaFusionData rnaData = new RnaFusionData(
                        items[COL_NAME], geneUp, geneDown, items[COL_CHR_UP], items[COL_CHR_DOWN],
                        Long.parseLong(items[COL_POS_UP]), Long.parseLong(items[COL_POS_DOWN]),
                        Byte.parseByte(items[COL_STRAND_UP]), Byte.parseByte(items[COL_STRAND_DOWN]),
                        Integer.parseInt(items[COL_JUNCT_RC]),Integer.parseInt(items[COL_SPAN_RC]), items[COL_SPLICE]);

                rnaDataList.add(rnaData);

                line = fileReader.readLine();
            }

        }
        catch(IOException e)
        {
            LOGGER.warn("failed to load sample RNA data file({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

    private String checkAlternateGeneName(final String geneName)
    {
        if(geneName.equals("AC005152.2"))
            return "SOX9-AS1";

        if(geneName.equals("AC016683.6"))
            return "PAX8-AS1";

        if(geneName.equals("AC007092.1"))
            return "LINC01122";

        if(geneName.toUpperCase().equals("C10ORF112"))
            return "MALRD1";

        if(geneName.equals("C5orf50"))
            return "SMIM23";

        if(geneName.equals("C10orf68"))
            return geneName.toUpperCase();

        if(geneName.equals("C17orf76-AS1"))
            return "FAM211A-AS1";

        if(geneName.equals("IGH@") || geneName.equals("IGH-@"))
            return "IGHJ6";

        if(geneName.equals("IGL@") || geneName.equals("IGL-@"))
            return "IGLC6";

        if(geneName.equals("MKLN1-AS1"))
            return "LINC-PINT";

        if(geneName.equals("PHF15"))
            return "JADE2";

        if(geneName.equals("PHF17"))
            return "JADE1";

        if(geneName.equals("RP11-134P9.1"))
            return "LINC01136";

        if(geneName.equals("RP11-973F15.1"))
            return "LINC01151";

        if(geneName.equals("RP11-115K3.2"))
            return "YWHAEP7";

        if(geneName.equals("RP11-3B12.1"))
            return "POT1-AS1";

        if(geneName.equals("RP11-199O14.1"))
            return "CASC20";

        if(geneName.equals("RP11-264F23.3"))
            return "CCND2-AS1";

        if(geneName.equals("RP11-93L9.1"))
            return "LINC01091";

        return geneName;
    }

    public void close()
    {
        closeBufferedWriter(mWriter);
    }

}
