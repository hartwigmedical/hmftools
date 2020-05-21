package com.hartwig.hmftools.linx.fusion.rna;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.EXON_RANK_MIN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWNSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UPSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.streamStr;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.linx.fusion.FusionFinder.checkFusionLogic;
import static com.hartwig.hmftools.linx.fusion.rna.RnaFusionAnnotator.findExonMatch;
import static com.hartwig.hmftools.linx.fusion.rna.RnaFusionAnnotator.isViableBreakend;
import static com.hartwig.hmftools.linx.fusion.rna.RnaFusionAnnotator.setReferenceFusionData;
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
import com.hartwig.hmftools.common.fusion.GeneAnnotation;
import com.hartwig.hmftools.common.fusion.GeneFusion;
import com.hartwig.hmftools.common.fusion.Transcript;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.fusion.FusionFinder;
import com.hartwig.hmftools.linx.fusion.FusionParameters;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvVarData;

public class RnaFusionMapper
{
    private String mSampleId;

    private final FusionFinder mFusionFinder;
    private final RnaMatchWriter mWriter;
    private final FusionParameters mFusionParams;
    private final EnsemblDataCache mGeneTransCache;
    private final Map<String, List<RnaFusionData>> mSampleRnaData;
    private final RnaFusionAnnotator mAnnotator;

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
        mAnnotator = new RnaFusionAnnotator(geneTransCache);
        mWriter = new RnaMatchWriter(outputDir);

        mFusionParams = new FusionParameters();
        mFusionParams.RequirePhaseMatch = false;
        mFusionParams.AllowExonSkipping = false;
    }

    public final Map<String, List<RnaFusionData>> getSampleRnaData() { return mSampleRnaData; }
    public final List<RnaFusionData> getSampleRnaData(final String sampleId) { return mSampleRnaData.get(sampleId); }
    public final RnaFusionAnnotator getAnnotator() { return mAnnotator; }

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
        final List<SvBreakend>[] viableBreakendPair = new List[] { Lists.newArrayList(), Lists.newArrayList() };
        final List<Transcript>[] viableTranscriptPair = new List[] { Lists.newArrayList(), Lists.newArrayList() };

        // transcripts on the correct side and orientation of the RNA boundary
        final List<Transcript>[] nearTranscriptPair = new List[] { Lists.newArrayList(), Lists.newArrayList() };
        final List<SvBreakend>[] nearBreakendPair = new List[] { Lists.newArrayList(), Lists.newArrayList() };

        // non-viable transcripts to be used if no others are found
        final List<Transcript>[] genicTranscriptPair = new List[] { Lists.newArrayList(), Lists.newArrayList() };
        final List<SvBreakend>[] genicBreakendPair = new List[] { Lists.newArrayList(), Lists.newArrayList() };

        boolean isExactRnaExon = rnaFusion.matchesKnownSpliceSite();

        for(int fs = FS_UPSTREAM; fs <= FS_DOWNSTREAM ; ++fs)
        {
            boolean isUpstream = (fs == 0);
            final String chromosome = rnaFusion.Chromosomes[fs];
            int rnaPosition = rnaFusion.Positions[fs];
            byte geneStrand = rnaFusion.Strands[fs];
            List<SvBreakend> viableBreakends = viableBreakendPair[fs];
            List<SvBreakend> nearBreakends = nearBreakendPair[fs];
            List<SvBreakend> genicBreakends = genicBreakendPair[fs];
            List<Transcript> viableTranscripts = viableTranscriptPair[fs];
            List<Transcript> nearTranscripts = nearTranscriptPair[fs];
            List<Transcript> genicTranscripts = genicTranscriptPair[fs];
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

                    if(correctLocation && mAnnotator.isTranscriptBreakendViableForRnaBoundary(
                            trans, isUpstream,  breakend.position(), rnaPosition, isExactRnaExon))
                    {
                        viableBreakends.add(breakend);
                        viableTranscripts.add(trans);
                    }
                }
            }
        }

        LNX_LOGGER.debug("rna fusion({}) breakend matches: upstream(viable={} near={} genic={}) downstream(viable={} near={} genic={})",
                rnaFusion.name(), viableBreakendPair[FS_UPSTREAM].size(), nearBreakendPair[FS_UPSTREAM].size(), genicBreakendPair[FS_UPSTREAM].size(),
                viableBreakendPair[FS_DOWNSTREAM].size(), nearBreakendPair[FS_DOWNSTREAM].size(), genicBreakendPair[FS_DOWNSTREAM].size());

        // run them through fusion logic (ie a pair of breakend lists), but don't require phase matching
        if(!viableBreakendPair[FS_UPSTREAM].isEmpty() && !viableBreakendPair[FS_DOWNSTREAM].isEmpty())
        {
            GeneFusion topCandidateFusion = null;
            SvBreakend topUpBreakend = null;
            SvBreakend topDownBreakend = null;
            boolean topCandidateFusionViable = false;

            for (int i = 0; i < viableBreakendPair[FS_UPSTREAM].size(); ++i)
            {
                final SvBreakend upBreakend = viableBreakendPair[FS_UPSTREAM].get(i);
                final Transcript upTrans = viableTranscriptPair[FS_UPSTREAM].get(i);

                if(upBreakend.getSV().isSglBreakend())
                    continue;

                for (int j = 0; j < viableBreakendPair[FS_DOWNSTREAM].size(); ++j)
                {
                    final SvBreakend downBreakend = viableBreakendPair[FS_DOWNSTREAM].get(j);
                    final Transcript downTrans = viableTranscriptPair[FS_DOWNSTREAM].get(j);

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
                if(!viableTranscriptPair[fs].isEmpty())
                {
                    isViable = true;
                    correctLocation = true;
                    transcriptList = viableTranscriptPair[fs];
                    breakendList = viableBreakendPair[fs];
                }
                else if(!nearTranscriptPair[fs].isEmpty())
                {
                    correctLocation = true;
                    transcriptList = nearTranscriptPair[fs];
                    breakendList = nearBreakendPair[fs];
                }
                else
                {
                    transcriptList = genicTranscriptPair[fs];
                    breakendList = genicBreakendPair[fs];
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

    private void setRnaFusionData(final RnaFusionData rnaFusion)
    {
        // find transcripts which match the RNA positions
        Map<String,RnaExonMatchData> transPhasesUp = Maps.newHashMap();
        Map<String,RnaExonMatchData> transPhasesDown = Maps.newHashMap();

        boolean isExactRnaExon = rnaFusion.matchesKnownSpliceSite();

        for(int fs = FS_UPSTREAM; fs <= 1; ++fs)
        {
            boolean isUpstream = (fs == 0);
            int rnaPosition = rnaFusion.Positions[fs];

            Map<String,RnaExonMatchData> transPhases = isUpstream ? transPhasesUp : transPhasesDown;

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

                    RnaExonMatchData exonMatchData = findExonMatch(transData, rnaPosition);

                    if(exonMatchData.BoundaryMatch)
                        rnaFusion.getExactMatchTransIds()[fs].add(transData.TransName);

                    if(exonMatchData.ExonFound)
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
            for (Map.Entry<String,RnaExonMatchData> entryUp : transPhasesUp.entrySet())
            {
                final String transIdUp = entryUp.getKey();
                final RnaExonMatchData exonDataUp = entryUp.getValue();

                if(isExactRnaExon && !rnaFusion.getExactMatchTransIds()[FS_UPSTREAM].contains(transIdUp))
                    continue;

                for (Map.Entry<String,RnaExonMatchData> entryDown : transPhasesDown.entrySet())
                {
                    final String transIdDown = entryDown.getKey();
                    final RnaExonMatchData exonDataDown = entryDown.getValue();

                    if(isExactRnaExon && !rnaFusion.getExactMatchTransIds()[FS_DOWNSTREAM].contains(transIdDown))
                        continue;

                    boolean phaseMatched = exonDataUp.ExonPhase == exonDataDown.ExonPhase;

                    if (phaseMatched && !rnaFusion.hasRnaPhasedFusion())
                    {
                        LNX_LOGGER.debug("rnaFusion({}) juncTypes(up={} down={}) transUp({}) transDown({} phase({}) matched({})",
                                rnaFusion.name(), rnaFusion.JunctionTypes[FS_UPSTREAM], rnaFusion.JunctionTypes[FS_DOWNSTREAM],
                                transIdUp, transIdDown, exonDataUp.ExonPhase, phaseMatched);

                        rnaFusion.setRnaPhasedFusionData(transIdUp, transIdDown);
                        rnaFusion.setExonData(FS_UPSTREAM, exonDataUp.ExonRank, exonDataUp.ExonPhase);
                        rnaFusion.setExonData(FS_DOWNSTREAM, exonDataDown.ExonRank, exonDataDown.ExonPhase);
                    }
                }
            }
        }

        setReferenceFusionData(mFusionFinder.getKnownFusionData(), rnaFusion);
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
