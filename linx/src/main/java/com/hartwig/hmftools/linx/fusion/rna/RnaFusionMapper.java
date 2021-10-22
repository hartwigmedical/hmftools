package com.hartwig.hmftools.linx.fusion.rna;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.EXON_RANK_MIN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.streamStr;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.MAX_UPSTREAM_DISTANCE_KNOWN;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.MAX_UPSTREAM_DISTANCE_OTHER;
import static com.hartwig.hmftools.linx.fusion.FusionFinder.checkFusionLogic;
import static com.hartwig.hmftools.linx.fusion.rna.RnaDataLoader.RNA_FUSION_SOURCE_ISOFOX;
import static com.hartwig.hmftools.linx.fusion.rna.RnaDataLoader.getRnaSourceDelimiter;
import static com.hartwig.hmftools.linx.fusion.rna.RnaDataLoader.loadRnaFusion;
import static com.hartwig.hmftools.linx.fusion.rna.RnaFusionAnnotator.checkRnaPhasedTranscripts;
import static com.hartwig.hmftools.linx.fusion.rna.RnaFusionAnnotator.findExonMatch;
import static com.hartwig.hmftools.linx.fusion.rna.RnaFusionAnnotator.isViableBreakend;
import static com.hartwig.hmftools.linx.fusion.rna.RnaFusionAnnotator.setReferenceFusionData;
import static com.hartwig.hmftools.linx.fusion.rna.RnaFusionAnnotator.unsplicedPositionMatch;
import static com.hartwig.hmftools.linx.fusion.rna.RnaJunctionType.KNOWN;
import static com.hartwig.hmftools.linx.fusion.rna.RnaJunctionType.NOT_SET;
import static com.hartwig.hmftools.linx.fusion.rna.RnaJunctionType.UNKNOWN;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.linx.gene.BreakendGeneData;
import com.hartwig.hmftools.linx.gene.BreakendTransData;
import com.hartwig.hmftools.common.utils.sv.StartEndPair;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.fusion.FusionFinder;
import com.hartwig.hmftools.linx.fusion.FusionParameters;
import com.hartwig.hmftools.linx.fusion.GeneFusion;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.commons.cli.CommandLine;

public class RnaFusionMapper
{
    private String mSampleId;

    private final FusionFinder mFusionFinder;
    private final RnaMatchWriter mWriter;
    private final FusionParameters mFusionParams;
    private final EnsemblDataCache mGeneTransCache;
    private final Map<String,List<RnaFusionData>> mSampleRnaData;
    private final RnaFusionAnnotator mAnnotator;

    private final List<GeneFusion> mDnaFusions;
    private final Map<GeneFusion,String> mDnaInvalidFusions;

    public static final String RNA_FUSIONS_FILE = "rna_fusions_file";
    public static final String RNA_FILE_SOURCE = "rna_file_source";

    public RnaFusionMapper(final String outputDir, final CommandLine cmdLineArgs, final EnsemblDataCache geneTransCache,
            FusionFinder fusionFinder, final List<GeneFusion> dnaFusions, final Map<GeneFusion,String> dnaInvalidFusions)
    {
        mSampleRnaData = Maps.newHashMap();
        mFusionFinder = fusionFinder;
        mGeneTransCache = geneTransCache;
        mDnaFusions = dnaFusions;
        mDnaInvalidFusions = dnaInvalidFusions;
        mAnnotator = new RnaFusionAnnotator(geneTransCache);

        final String fileSource = cmdLineArgs.getOptionValue(RNA_FILE_SOURCE, RNA_FUSION_SOURCE_ISOFOX);

        if(cmdLineArgs != null)
        {
            final String rnaDataFile = cmdLineArgs.getOptionValue(RNA_FUSIONS_FILE);
            loadSampleRnaData(fileSource, rnaDataFile);
        }

        mWriter = new RnaMatchWriter(outputDir, fileSource);

        mFusionParams = new FusionParameters();
        mFusionParams.RequirePhaseMatch = false;
        mFusionParams.AllowExonSkipping = false;
    }

    public final Map<String, List<RnaFusionData>> getSampleRnaData() { return mSampleRnaData; }

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

            findMatchingBreakendTranscripts(rnaFusion, chrBreakendMap);

            setDnaFusionMatch(rnaFusion);

            setReferenceFusionData(mFusionFinder.getKnownFusionCache(), rnaFusion);

            rnaFusion.setFusionClusterChainInfo();

            mWriter.writeRnaMatchData(rnaFusion);
        }

        // move from consideration to de-link RNA data from SV types
        mSampleRnaData.remove(mSampleId);
    }

    private void findMatchingBreakendTranscripts(final RnaFusionData rnaFusion, Map<String, List<SvBreakend>> chrBreakendMap)
    {
        /* Matching and annotation logic:
            - find all breakends in the RNA up and down gene
            - for them, find the any transcripts which a) have the exon boundary in the RNA position AND
            - b) are in the correct relative position:
                - upstream: at or after the RNA boundary down to the start of the next exon
                - downstream: at or before the RNA boundary up to the start of the preceding exon
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
        final StartEndPair<List<SvBreakend>> viableBreakendPair = new StartEndPair<>(Lists.newArrayList(), Lists.newArrayList());
        final StartEndPair<List<BreakendTransData>> viableTranscriptPair = new StartEndPair<>(Lists.newArrayList(), Lists.newArrayList());

        // transcripts on the correct side and orientation of the RNA boundary
        final StartEndPair<List<BreakendTransData>> nearTranscriptPair = new StartEndPair<>(Lists.newArrayList(), Lists.newArrayList());
        final StartEndPair<List<SvBreakend>> nearBreakendPair = new StartEndPair<>(Lists.newArrayList(), Lists.newArrayList());

        boolean requireExactMatch = rnaFusion.JunctionTypes[FS_UP] == UNKNOWN || rnaFusion.JunctionTypes[FS_DOWN] == UNKNOWN;

        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            boolean isUpstream = (fs == 0);
            final String chromosome = rnaFusion.Chromosomes[fs];
            int rnaPosition = rnaFusion.Positions[fs];
            byte rnaOrient = rnaFusion.Orientations[fs];
            byte requiredGeneStrand = isUpstream ? rnaOrient : (byte)-rnaOrient;
            String geneName = rnaFusion.GeneNames[fs];
            boolean isKnownJunction = rnaFusion.JunctionTypes[fs] == KNOWN;
            int maxPreGeneDistance = isKnownJunction ? MAX_UPSTREAM_DISTANCE_KNOWN : MAX_UPSTREAM_DISTANCE_OTHER;

            final List<SvBreakend> breakendList = chrBreakendMap.get(chromosome);

            if(breakendList == null)
                continue;

            for(final SvBreakend breakend : breakendList)
            {
                if(breakend.orientation() != rnaOrient)
                     continue;

                final SvVarData var = breakend.getSV();

                if(var.isInferredSgl())
                    continue;

                final List<BreakendGeneData> svGenesList = var.getGenesList(breakend.usesStart());

                if(svGenesList.isEmpty())
                    continue;

                // breakend must be within a max distance of the RNA breakend or if it's unspliced, at its exact location
                if(requireExactMatch)
                {
                    if(!unsplicedPositionMatch(breakend, rnaPosition))
                        continue;
                }
                else
                {
                    if(!isViableBreakend(breakend, rnaPosition, requiredGeneStrand, isUpstream))
                        continue;
                }

                final List<BreakendGeneData> genesList = Lists.newArrayList();

                for(BreakendGeneData gene : svGenesList)
                {
                    if(!geneName.isEmpty() && !geneName.equals(gene.geneName()))
                        continue;

                    if(!requireExactMatch)
                    {
                        // breakends cannot be too far upstream from any gene being considered
                        if(gene.strand() == POS_STRAND && gene.GeneData.GeneStart - breakend.position() > maxPreGeneDistance)
                            continue;
                        else if(gene.strand() == NEG_STRAND && breakend.position() - gene.GeneData.GeneEnd > maxPreGeneDistance)
                            continue;
                    }

                    genesList.add(gene);
                }

                if(genesList.isEmpty())
                    continue;

                // check that breakend has correct orientation and position relative to RNA breakend
                // boolean correctLocation = isViableBreakend(breakend, rnaPosition, requiredGeneStrand, isUpstream);

                // check whether any of the breakend's transcripts match the exon (exact or nearest) of the RNA fusion breakpoint
                for(final BreakendGeneData gene : genesList)
                {
                    for(final BreakendTransData trans : gene.transcripts())
                    {
                        if(isKnownJunction && !rnaFusion.getExactMatchTransIds(fs).contains(trans.transName()))
                            continue;

                        nearBreakendPair.get(fs).add(breakend);
                        nearTranscriptPair.get(fs).add(trans);

                        if(mAnnotator.isTranscriptBreakendViableForRnaBoundary(
                                trans, isUpstream, breakend.position(), rnaPosition, isKnownJunction))
                        {
                            viableBreakendPair.get(fs).add(breakend);
                            viableTranscriptPair.get(fs).add(trans);
                        }
                    }
                }
            }
        }

        LNX_LOGGER.debug("rna fusion({}) breakend matches: upstream(viable={} near={}) downstream(viable={} near={})",
                rnaFusion.name(), viableBreakendPair.get(FS_UP).size(), nearBreakendPair.get(FS_UP).size(),
                viableBreakendPair.get(FS_DOWN).size(), nearBreakendPair.get(FS_DOWN).size());

        // run them through fusion logic (ie a pair of breakend lists), but don't require phase matching
        if(!viableBreakendPair.get(FS_UP).isEmpty() && !viableBreakendPair.get(FS_DOWN).isEmpty())
        {
            GeneFusion topCandidateFusion = null;
            SvBreakend topUpBreakend = null;
            SvBreakend topDownBreakend = null;
            boolean topCandidateFusionViable = false;

            for (int i = 0; i < viableBreakendPair.get(FS_UP).size(); ++i)
            {
                final SvBreakend upBreakend = viableBreakendPair.get(FS_UP).get(i);
                final BreakendTransData upTrans = viableTranscriptPair.get(FS_UP).get(i);

                if(upBreakend.getSV().isSglBreakend())
                    continue;

                for (int j = 0; j < viableBreakendPair.get(FS_DOWN).size(); ++j)
                {
                    final SvBreakend downBreakend = viableBreakendPair.get(FS_DOWN).get(j);
                    final BreakendTransData downTrans = viableTranscriptPair.get(FS_DOWN).get(j);

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
                        FS_UP, topCandidateFusion.upstreamTrans(), topUpBreakend,
                        true, true,  0);

                rnaFusion.setTranscriptData(
                        FS_DOWN, topCandidateFusion.downstreamTrans(), topDownBreakend,
                        true, true,0);

                rnaFusion.setViableFusion(topCandidateFusionViable, topCandidateFusion.phaseMatched());
            }
        }
        else
        {
            // select the closest breakend's transcript
            for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
            {
                boolean isUpstream = (fs == 0);
                int rnaPosition = rnaFusion.Positions[fs];

                List<BreakendTransData> transcriptList;
                List<SvBreakend> breakendList;
                boolean isViable = false;
                boolean correctLocation = false;

                // use the viable transcripts if present, otherwise the nearest
                if(!viableTranscriptPair.get(fs).isEmpty())
                {
                    isViable = true;
                    correctLocation = true;
                    transcriptList = viableTranscriptPair.get(fs);
                    breakendList = viableBreakendPair.get(fs);
                }
                else if(!nearTranscriptPair.get(fs).isEmpty())
                {
                    correctLocation = true;
                    transcriptList = nearTranscriptPair.get(fs);
                    breakendList = nearBreakendPair.get(fs);
                }
                else
                {
                    continue;
                    // transcriptList = genicTranscriptPair.get(fs);
                    // breakendList = genicBreakendPair.get(fs);
                }

                BreakendTransData closestTrans = null;
                SvBreakend closestBreakend = null;
                int closestDistance = 0;

                for (int j = 0; j < transcriptList.size(); ++j)
                {
                    final BreakendTransData trans = transcriptList.get(j);
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
                        String geneId = closestTrans.gene().geneId();
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
    }

    private void setDnaFusionMatch(final RnaFusionData rnaFusion)
    {
        DnaRnaMatchType matchType = DnaRnaMatchType.NONE;
        boolean reportableFusion = false;

        boolean isUnspliced = rnaFusion.JunctionTypes[FS_UP] == UNKNOWN && rnaFusion.JunctionTypes[FS_DOWN] == UNKNOWN;

        final BreakendTransData transUp = rnaFusion.getMatchedTranscripts()[FS_UP];
        final BreakendTransData transDown = rnaFusion.getMatchedTranscripts()[FS_DOWN];

        for(final GeneFusion dnaFusion : mDnaFusions)
        {
            if(transUp != null && dnaFusion.upstreamTrans().gene().id() == transUp.gene().id()
            && transDown != null && dnaFusion.downstreamTrans().gene().id() == transDown.gene().id())
            {
                matchType = DnaRnaMatchType.SVS;
                reportableFusion = dnaFusion.reportable();

                if(!rnaFusion.GeneNames[FS_UP].isEmpty() && !rnaFusion.GeneNames[FS_DOWN].isEmpty()
                && (!dnaFusion.upstreamTrans().geneName().equals(rnaFusion.GeneNames[FS_UP])
                || !dnaFusion.downstreamTrans().geneName().equals(rnaFusion.GeneNames[FS_DOWN])))
                {
                    LNX_LOGGER.debug("genePair rna({}-{}) differs from dna({}-{}) for same SVs({} & {})",
                            rnaFusion.GeneNames[FS_UP], rnaFusion.GeneNames[FS_DOWN],
                            dnaFusion.upstreamTrans().geneName(), dnaFusion.downstreamTrans().geneName(),
                            transUp.gene().id(), transDown.gene().id());

                    // override the gene selection for downstream analysis
                    rnaFusion.GeneNames[FS_UP] = dnaFusion.upstreamTrans().geneName();
                    rnaFusion.GeneNames[FS_DOWN] = dnaFusion.downstreamTrans().geneName();
                }

                break;
            }

            if(!dnaFusion.upstreamTrans().geneName().equals(rnaFusion.GeneNames[FS_UP])
            || !dnaFusion.downstreamTrans().geneName().equals(rnaFusion.GeneNames[FS_DOWN]))
            {
                continue;
            }

            matchType = DnaRnaMatchType.GENES;
            reportableFusion = dnaFusion.reportable();
        }

        if(matchType != DnaRnaMatchType.NONE)
        {
            rnaFusion.setDnaFusionMatch(matchType, "", reportableFusion);
            return;
        }

        if(isUnspliced)
            return; // cannot match an invalid chained fusion

        // search amongst the invalid chained fusions for a match on SVs
        for(Map.Entry<GeneFusion,String> entry : mDnaInvalidFusions.entrySet())
        {
            final GeneFusion fusion = entry.getKey();

            if(transUp != null && fusion.upstreamTrans().gene().id() == transUp.gene().id()
            && transDown != null && fusion.downstreamTrans().gene().id() == transDown.gene().id())
            {
                rnaFusion.setDnaFusionMatch(DnaRnaMatchType.INVALID, entry.getValue(), false);
                return;
            }
        }

        // then by gene pair
        for(Map.Entry<GeneFusion,String> entry : mDnaInvalidFusions.entrySet())
        {
            final GeneFusion fusion = entry.getKey();

            if(fusion.name().equals(rnaFusion.name()))
            {
                rnaFusion.setDnaFusionMatch(DnaRnaMatchType.INVALID, entry.getValue(), false);
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
        double currentPosDiff = (abs(rnaFusion.Positions[FS_UP] - beCurrentStart.position())
                + abs(rnaFusion.Positions[FS_DOWN] - beCurrentEnd.position())) * 0.5;

        double candidatePosDiff = (abs(rnaFusion.Positions[FS_UP] - beCandidateStart.position())
                + abs(rnaFusion.Positions[FS_DOWN] - beCandidateEnd.position())) * 0.5;

        return candidatePosDiff < currentPosDiff;
    }

    private void setRnaFusionData(final RnaFusionData rnaFusion)
    {
        // correct gene names
        mAnnotator.correctGeneNames(mFusionFinder.getKnownFusionCache(), rnaFusion);

        // find transcripts which match the RNA positions
        for(int fs = FS_UP; fs <= 1; ++fs)
        {
            int rnaPosition = rnaFusion.Positions[fs];

            if(rnaFusion.GeneIds[fs].isEmpty())
                continue;

            GeneData geneData = mGeneTransCache.getGeneDataById(rnaFusion.GeneIds[fs]);

            if(geneData == null)
            {
                LNX_LOGGER.warn("sample({}) rnaFusion({}) {} gene not found", mSampleId, rnaFusion.name(), streamStr(fs));
                continue;
            }

            rnaFusion.Strands[fs] = geneData.Strand;

            List<TranscriptData> transDataList = mGeneTransCache.getTranscripts(geneData.GeneId);

            if(transDataList != null)
            {
                for(final TranscriptData transData : transDataList)
                {
                    RnaExonMatchData exonMatchData = findExonMatch(transData, rnaPosition);

                    if(exonMatchData.ExonFound)
                    {
                        rnaFusion.getTransExonData(fs).add(exonMatchData);

                        // override for external fusion source data which doesn't set this or set it correctly
                        if(exonMatchData.BoundaryMatch && rnaFusion.JunctionTypes[fs] == NOT_SET)
                            rnaFusion.JunctionTypes[fs] = KNOWN;
                    }
                }
            }

            if(rnaFusion.JunctionTypes[fs] == NOT_SET)
                rnaFusion.JunctionTypes[fs] = UNKNOWN;

            LNX_LOGGER.debug("rnaFusion({}) juncType({}) {} position({}) matched {} transcript exons",
                    rnaFusion.name(), rnaFusion.JunctionTypes[fs], streamStr(fs), rnaPosition, rnaFusion.getTransExonData(fs).size());
        }

        checkRnaPhasedTranscripts(rnaFusion);
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
                RnaFusionData rnaData = loadRnaFusion(source, recordCount, data, fieldIndexMap);
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
