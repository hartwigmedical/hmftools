package com.hartwig.hmftools.linx.fusion;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.EXON_DEL_DUP;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_PROMISCUOUS;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.NONE;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_3;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_5;

import static com.hartwig.hmftools.common.fusion.KnownFusionCache.KNOWN_FUSIONS_FILE;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.ENHANCER;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.tickPhaseForward;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.REQUIRED_BIOTYPES;
import static com.hartwig.hmftools.linx.fusion.FusionReportability.checkProteinDomains;
import static com.hartwig.hmftools.linx.fusion.FusionReportability.determineReportability;
import static com.hartwig.hmftools.linx.fusion.FusionReportability.findTopPriorityFusion;
import static com.hartwig.hmftools.linx.fusion.FusionReportability.validProteinDomains;
import static com.hartwig.hmftools.linx.fusion.ReportableReason.PROTEIN_DOMAINS;
import static com.hartwig.hmftools.linx.fusion.ReportableReason.OK;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.linx.gene.BreakendGeneData;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.linx.gene.BreakendTransData;
import com.hartwig.hmftools.common.gene.TranscriptProteinData;
import com.hartwig.hmftools.common.gene.TranscriptRegionType;

import org.apache.commons.cli.Options;

public class FusionFinder
{
    private final KnownFusionCache mKnownFusionCache;

    private final EnsemblDataCache mGeneTransCache;
    public int mNextFusionId;

    private static boolean mLogInvalidReasons;

    public FusionFinder(final EnsemblDataCache geneTransCache, final KnownFusionCache knownFusionCache)
    {
        mGeneTransCache = geneTransCache;
        mKnownFusionCache = knownFusionCache;

        mNextFusionId = 0;

        FusionReportability.populateRequiredProteins();

        mLogInvalidReasons = false;
    }

    public boolean hasValidConfigData() { return mKnownFusionCache.hasValidData(); }
    public void setLogInvalidReasons(boolean toggle) { mLogInvalidReasons = toggle; }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(KNOWN_FUSIONS_FILE, true, "Known fusion file");
    }

    public final KnownFusionCache getKnownFusionCache() { return mKnownFusionCache; }

    public static final String INVALID_REASON_ORIENTATION = "Orientation";
    public static final String INVALID_REASON_PHASING = "Unphased";
    public static final String INVALID_REASON_CODING_TYPE = "Coding";

    public void reset() { mNextFusionId = 0; }

    public final List<GeneFusion> findFusions(
            final List<BreakendGeneData> breakendGenes1, final List<BreakendGeneData> breakendGenes2, final FusionParameters params)
    {
        final List<GeneFusion> potentialFusions = Lists.newArrayList();

        if(!hasValidConfigData())
            return potentialFusions;

        for (final BreakendGeneData startGene : breakendGenes1)
        {
            // left is upstream, right is downstream
            boolean startUpstream = startGene.isUpstream();
            boolean startIsIgRegion = mKnownFusionCache.withinIgRegion(startGene.chromosome(), startGene.position());

            for (final BreakendGeneData endGene : breakendGenes2)
            {
                boolean endUpstream = endGene.isUpstream();
                boolean endIsIgRegion = mKnownFusionCache.withinIgRegion(endGene.chromosome(), endGene.position());

                if(startIsIgRegion || endIsIgRegion)
                {
                    if(startIsIgRegion && endIsIgRegion)
                        continue;

                    checkIgFusion(startGene, endGene, potentialFusions);
                    continue;
                }

                if (startUpstream == endUpstream)
                {
                    if(params.InvalidReasons != null && !params.InvalidReasons.contains(INVALID_REASON_ORIENTATION))
                        params.InvalidReasons.add(INVALID_REASON_ORIENTATION);

                    continue;
                }

                final BreakendGeneData upGene = startUpstream ? startGene : endGene;
                final BreakendGeneData downGene = !startUpstream ? startGene : endGene;

                boolean knownPair = mKnownFusionCache.hasKnownFusion(upGene.GeneName, downGene.GeneName);
                boolean knownUnmappable3Pair = mKnownFusionCache.hasKnownUnmappable3Fusion(upGene.GeneName, downGene.GeneName);

                for(final BreakendTransData upstreamTrans : upGene.transcripts())
                {
                    if(!isValidUpstreamTranscript(upstreamTrans, !knownPair, params.RequireUpstreamBiotypes))
                    {
                        logInvalidReasonInfo(upstreamTrans, null, INVALID_REASON_CODING_TYPE, "invalid up trans");
                        continue;
                    }

                    for(final BreakendTransData downstreamTrans : downGene.transcripts())
                    {
                        GeneFusion fusion;

                        if(knownUnmappable3Pair)
                        {
                            // no further conditions
                            fusion = new GeneFusion(upstreamTrans, downstreamTrans, true);
                            fusion.setId(mNextFusionId++);
                        }
                        else
                        {
                            if(!isValidDownstreamTranscript(downstreamTrans))
                            {
                                logInvalidReasonInfo(upstreamTrans, downstreamTrans, INVALID_REASON_CODING_TYPE, "invalid down trans");
                                continue;
                            }

                            boolean exonDelDupCandidate = upstreamTrans.transName().equals(downstreamTrans.transName())
                                    && mKnownFusionCache.isExonDelDupTrans(upstreamTrans.transName());

                            fusion = transcriptPairHasFusion(upstreamTrans, downstreamTrans, params, knownPair, exonDelDupCandidate);
                        }

                        if(fusion == null)
                            continue;

                        fusion.setId(mNextFusionId++);
                        setKnownFusionType(fusion);
                        potentialFusions.add(fusion);
                    }
                }
            }
        }

        return potentialFusions;
    }

    private static void logInvalidReasonInfo(final BreakendTransData trans1, final BreakendTransData trans2, final String reasonType, final String reason)
    {
        if(!mLogInvalidReasons)
            return;

        if(trans2 == null)
            LNX_LOGGER.trace("transcript({}:{}) invalid({}: {})", trans1.geneName(), trans1.transName(), reasonType, reason);
        else
            LNX_LOGGER.trace("transcripts({}:{} and {}:{}) invalid({}: {})",
                    trans1.geneName(), trans1.transName(), trans2.geneName(), trans2.transName(), reasonType, reason);
    }

    public static boolean validFusionTranscript(final BreakendTransData transcript)
    {
        if(transcript.isUpstream())
            return isValidUpstreamTranscript(transcript, true, false);
        else
            return isValidDownstreamTranscript(transcript);
    }

    private static boolean isValidUpstreamTranscript(
            final BreakendTransData transcript, boolean requireUpstreamDisruptive, boolean requireUpstreamBiotypes)
    {
        // check any conditions which would preclude this transcript being a part of a fusion no matter the other end
        if(transcript.isPromoter())
            return false;

        if(requireUpstreamDisruptive && !transcript.isDisruptive())
            return false;

        if(requireUpstreamBiotypes && !REQUIRED_BIOTYPES.contains(transcript.bioType()))
            return false;

        return true;
    }

    private static boolean isValidDownstreamTranscript(final BreakendTransData transcript)
    {
        if(transcript.postCoding())
            return false;

        if(transcript.nonCoding())
            return false;

        if(transcript.exonCount() == 1)
            return false;

        return true;
    }

    public static GeneFusion checkFusionLogic(final BreakendTransData upstreamTrans, final BreakendTransData downstreamTrans, final FusionParameters params)
    {
        if(!isValidUpstreamTranscript(upstreamTrans, true, params.RequireUpstreamBiotypes))
        {
            logInvalidReasonInfo(upstreamTrans, downstreamTrans, INVALID_REASON_CODING_TYPE, "invalid up trans");
            return null;
        }
        else if(!isValidDownstreamTranscript(downstreamTrans))
        {
            logInvalidReasonInfo(upstreamTrans, downstreamTrans, INVALID_REASON_CODING_TYPE, "invalid down trans");
            return null;
        }

        return transcriptPairHasFusion(upstreamTrans, downstreamTrans, params, false, false);
    }

    private static GeneFusion transcriptPairHasFusion(
            final BreakendTransData upstreamTrans, final BreakendTransData downstreamTrans, final FusionParameters params,
            boolean isKnownPair, boolean exonDelDupCandidate)
    {
        // see SV Fusions document for permitted combinations
        boolean checkExactMatch = false;

        if(upstreamTrans.preCoding())
        {
            if(upstreamTrans.isExonic() && !downstreamTrans.isExonic())
            {
                logInvalidReasonInfo(upstreamTrans, downstreamTrans, INVALID_REASON_CODING_TYPE, "precoding exonic to non-exonic");
                return null;
            }
            else if(downstreamTrans.isCoding() && !isKnownPair)
            {
                logInvalidReasonInfo(upstreamTrans, downstreamTrans, INVALID_REASON_CODING_TYPE, "pre-coding to coding");
                return null;
            }
            else if(downstreamTrans.preCoding() && upstreamTrans.gene().StableId.equals(downstreamTrans.gene().StableId))
            {
                // skip pre-coding to pre-coding within the same gene
                logInvalidReasonInfo(upstreamTrans, downstreamTrans, INVALID_REASON_CODING_TYPE, "pre-coding to pre-coding");
                return null;
            }
        }
        else if(upstreamTrans.isCoding())
        {
            if(downstreamTrans.nonCoding())
            {
                logInvalidReasonInfo(upstreamTrans, downstreamTrans, INVALID_REASON_CODING_TYPE, "coding to non-coding");
                return null;
            }

            if(downstreamTrans.preCoding() && !isKnownPair)
            {
                logInvalidReasonInfo(upstreamTrans, downstreamTrans, INVALID_REASON_CODING_TYPE, "coding to pre-coding");
                return null;
            }

            if(upstreamTrans.isExonic())
            {
                if(!downstreamTrans.isExonic() && (!params.AllowExonSkipping || !downstreamTrans.isIntronic()))
                {
                    logInvalidReasonInfo(upstreamTrans, downstreamTrans, INVALID_REASON_CODING_TYPE, "coding exonic to non-exonic");
                    return null;
                }

                if(upstreamTrans.gene().id() != downstreamTrans.gene().id())
                {
                    logInvalidReasonInfo(upstreamTrans, downstreamTrans, INVALID_REASON_CODING_TYPE, "up coding exonic diff SVs");
                    return null;
                }

                // coding exon to coding exon will require phase adjustments to be exact
                checkExactMatch = downstreamTrans.isExonic();
            }
        }
        else if(upstreamTrans.nonCoding())
        {
            if(upstreamTrans.isExonic() && !downstreamTrans.isExonic() && !isKnownPair)
            {
                logInvalidReasonInfo(upstreamTrans, downstreamTrans, INVALID_REASON_CODING_TYPE, "up non-coding exonic to down non-exonic");
                return null;
            }
            else if(downstreamTrans.isCoding() && !isKnownPair)
            {
                logInvalidReasonInfo(upstreamTrans, downstreamTrans, INVALID_REASON_CODING_TYPE, "up non-coding to down-coding");
                return null;
            }
        }

        if (isIrrelevantSameGene(upstreamTrans, downstreamTrans))
        {
            logInvalidReasonInfo(upstreamTrans, downstreamTrans, INVALID_REASON_CODING_TYPE, "irrelevant fusion");
            return null;
        }

        boolean phaseMatched = false;
        int phaseExonsSkippedUp = 0;
        int phaseExonsSkippedDown = 0;

        if(!checkExactMatch)
        {
            // all fusions to downstream exons may be excluded, but for now definitely exclude those which end in the last exon
            if(downstreamTrans.isExonic() && downstreamTrans.ExonDownstream == downstreamTrans.exonCount() && !downstreamTrans.preCoding())
            {
                logInvalidReasonInfo(upstreamTrans, downstreamTrans, INVALID_REASON_CODING_TYPE, "downstream last exon");
                return null;
            }
        }

        if(checkExactMatch)
        {
            phaseMatched = exonToExonInPhase(upstreamTrans, downstreamTrans);

            if(!phaseMatched)
            {
                logInvalidReasonInfo(upstreamTrans, downstreamTrans, INVALID_REASON_PHASING, "exon-exon inexact phasing");
            }

            if(phaseMatched || !params.RequirePhaseMatch)
            {
                return new GeneFusion(upstreamTrans, downstreamTrans, phaseMatched);
            }
        }
        else
        {
            // just check for a phasing match
            if(upstreamTrans.isExonic() && downstreamTrans.isExonic())
            {
                phaseMatched = ((upstreamTrans.preCoding() || upstreamTrans.nonCoding()) && downstreamTrans.preCoding())
                        || (upstreamTrans.postCoding() && downstreamTrans.postCoding());
            }
            else
            {
                phaseMatched = upstreamTrans.Phase == downstreamTrans.Phase;
            }

            if(!phaseMatched && params.AllowExonSkipping
            && (!upstreamTrans.gene().StableId.equals(downstreamTrans.gene().StableId) || exonDelDupCandidate))
            {
                // check for a match within the alternative phasings from upstream and downstream of the breakend
                for (Map.Entry<Integer, Integer> altPhasing : upstreamTrans.getAlternativePhasing().entrySet())
                {
                    if (altPhasing.getKey() == downstreamTrans.Phase)
                    {
                        phaseMatched = true;
                        phaseExonsSkippedUp = altPhasing.getValue();
                        break;
                    }
                }

                if(!phaseMatched)
                {
                    for (Map.Entry<Integer, Integer> altPhasing : downstreamTrans.getAlternativePhasing().entrySet())
                    {
                        if (altPhasing.getKey() == upstreamTrans.Phase)
                        {
                            phaseMatched = true;
                            phaseExonsSkippedDown = altPhasing.getValue();
                            break;
                        }
                    }

                    if(!phaseMatched)
                    {
                        for (Map.Entry<Integer, Integer> altPhasingUp : upstreamTrans.getAlternativePhasing().entrySet())
                        {
                            for(Map.Entry<Integer, Integer> altPhasingDown : downstreamTrans.getAlternativePhasing().entrySet())
                            {
                                if(altPhasingUp.getKey() == altPhasingDown.getKey())
                                {
                                    phaseMatched = true;
                                    phaseExonsSkippedUp = altPhasingUp.getValue();
                                    phaseExonsSkippedDown = altPhasingDown.getValue();
                                    break;
                                }
                            }
                        }
                    }
                }

                // look for any match between the 2 sets of alt-phasings
            }

            // moving past an exon into the intron to phase-match is not considered as as skipping an exon
            if(upstreamTrans.isExonic() && phaseExonsSkippedUp == 1)
                phaseExonsSkippedUp = 0;

            if(downstreamTrans.isExonic() && phaseExonsSkippedDown == 1)
                phaseExonsSkippedDown = 0;

            if(!phaseMatched)
            {
                logInvalidReasonInfo(upstreamTrans, downstreamTrans, INVALID_REASON_PHASING, "inexact unphased");
            }

            if(phaseMatched || !params.RequirePhaseMatch)
            {
                GeneFusion fusion = new GeneFusion(upstreamTrans, downstreamTrans, phaseMatched);
                fusion.setExonsSkipped(phaseExonsSkippedUp, phaseExonsSkippedDown);
                return fusion;
            }
        }

        return null;
    }

    private static boolean exonToExonInPhase(final BreakendTransData upTrans, final BreakendTransData downTrans)
    {
        // check phasing and offset since exon start or coding start
        int upPhase = upTrans.ExonicBasePhase;
        int downPhase = downTrans.ExonicBasePhase;

        if(upPhase == PHASE_NONE && downPhase == PHASE_NONE)
            return true;

        return tickPhaseForward(upPhase) == downPhase;
    }

    private void checkIgFusion(final BreakendGeneData startGene, final BreakendGeneData endGene, final List<GeneFusion> potentialFusions)
    {
        // These are allowed to fuse with 5’UTR splice donors from 100kb upstream. Phasing is assumed to be -1
        // In some specific cases, allow fusions in the 3’UTR region or even downstream of the gene (eg BCL2 (common in Folicular Lymphomas[PP1])
        // in which case the 3' gene's orientation must be facing back upstream into the coding region

        // orientation in the IG region is no longer checked
        boolean startIsIgGene = mKnownFusionCache.withinIgRegion(startGene.chromosome(), startGene.position()); // matchesIgGene startGene.orientation()
        boolean endIsIgGene = !startIsIgGene && mKnownFusionCache.withinIgRegion(endGene.chromosome(), endGene.position()); // endGene.orientation()

        if(!startIsIgGene && !endIsIgGene)
            return;

        final BreakendGeneData igGene = startIsIgGene ? startGene : endGene;
        final BreakendGeneData downGene = startIsIgGene ? endGene : startGene;

        KnownFusionType knownType = NONE;

        final List<BreakendTransData> candidateTranscripts = downGene.transcripts().stream()
                .filter(x -> x.CodingBases == 0 && !x.isUpstream())
                .collect(Collectors.toList());

        KnownFusionData kfData = mKnownFusionCache.getDataByType(IG_KNOWN_PAIR).stream()
                .filter(x -> x.ThreeGene.equals(downGene.GeneName))
                .filter(x -> x.withinIgRegion(igGene.chromosome(), igGene.position()))
                .findFirst().orElse(null);

        if(kfData != null)
        {
            // a known IG-partner gene
            if(kfData.downstreamDistance(FS_DOWN) > 0)
            {
                // must face back to the coding region if 3' UTR or downstream of the gene
                candidateTranscripts.addAll(downGene.transcripts().stream()
                        .filter(x -> x.codingType().equals(UTR_3P) && x.isUpstream())
                        .collect(Collectors.toList()));
            }

            knownType = IG_KNOWN_PAIR;
        }
        else
        {
            kfData = mKnownFusionCache.getDataByType(IG_PROMISCUOUS).stream()
                    .filter(x -> x.withinIgRegion(igGene.chromosome(), igGene.position()))
                    .findFirst().orElse(null);

            // check within the promiscuous region bounds
            if(kfData == null)
                return;

            knownType = IG_PROMISCUOUS;
        }

        final BreakendTransData upTrans = generateIgTranscript(igGene, kfData);

        if(!candidateTranscripts.isEmpty())
        {
            for(final BreakendTransData downTrans : candidateTranscripts)
            {
                GeneFusion fusion = new GeneFusion(upTrans, downTrans, true);
                fusion.setId(mNextFusionId++);
                fusion.setKnownType(knownType);
                potentialFusions.add(fusion);
            }
        }
    }

    private BreakendTransData generateIgTranscript(final BreakendGeneData gene, final KnownFusionData knownFusionData)
    {
        TranscriptData transData = new TranscriptData(
                0, String.format("@%s", knownFusionData.FiveGene), gene.StableId, false, gene.Strand,
                knownFusionData.igRegion().start(), knownFusionData.igRegion().end(), null, null, "");

        BreakendTransData transcript = new BreakendTransData(
                gene, transData,  -1, -1, PHASE_NONE, PHASE_NONE, 0, 0);

        transcript.setCodingType(ENHANCER);
        transcript.setRegionType(TranscriptRegionType.IG);
        transcript.setIsDisruptive(true);
        return transcript;
    }

    public static boolean isIrrelevantSameGene(final BreakendTransData upTrans, final BreakendTransData downTrans)
    {
        if(!upTrans.geneName().equals(downTrans.geneName()))
            return false;

        // skip fusions between different transcripts in the same gene,
        if (!upTrans.transName().equals(downTrans.transName()))
            return true;

        if(upTrans.nonCoding())
            return true;

        // skip fusions within the same intron
        if(upTrans.isIntronic() && downTrans.isIntronic() && upTrans.ExonUpstream == downTrans.ExonUpstream)
            return true;

        return false;
    }


    public GeneFusion findTopReportableFusion(final List<GeneFusion> fusions)
    {
        // find candidate reportable fusions, and then set protein domain info for these
        final List<GeneFusion> candidateReportable = Lists.newArrayList();

        for(GeneFusion fusion : fusions)
        {
            ReportableReason reason = determineReportability(fusion);
            fusion.setReportableReason(reason);

            if(reason != OK)
                continue;

            candidateReportable.add(fusion);
        }

        if(candidateReportable.isEmpty())
        {
            return findTopPriorityFusion(fusions);
        }

        GeneFusion reportableFusion = findTopPriorityFusion(candidateReportable);

        if(reportableFusion != null)
        {
            if(checkProteinDomains(reportableFusion.knownType()))
            {
                // check impact on protein regions
                setFusionProteinFeatures(reportableFusion);

                if(validProteinDomains(reportableFusion))
                    reportableFusion.setReportable(true);
                else
                    reportableFusion.setReportableReason(PROTEIN_DOMAINS);
            }
            else
            {
                reportableFusion.setReportable(true);
            }
        }

        return reportableFusion;
    }

    public void setFusionProteinFeatures(GeneFusion fusion)
    {
        if(fusion.proteinFeaturesSet())
            return;

        fusion.setProteinFeaturesSet();

        final BreakendTransData downTrans = fusion.downstreamTrans();

        if(downTrans.nonCoding())
            return;

        final List<TranscriptProteinData> transProteinData = mGeneTransCache.getTranscriptProteinDataMap().get(downTrans.transId());

        if(transProteinData == null || transProteinData.isEmpty())
            return;

        List<String> processedFeatures = Lists.newArrayList();

        for(int i = 0; i < transProteinData.size(); ++i)
        {
            final TranscriptProteinData pfData = transProteinData.get(i);
            final String feature = pfData.HitDescription;

            if(processedFeatures.contains(feature))
                continue;

            // find start and end across all entries matching this feature
            int featureStart = pfData.SeqStart;
            int featureEnd = pfData.SeqEnd;

            for(int j = i+1; j < transProteinData.size(); ++j)
            {
                if(transProteinData.get(j).HitDescription.equals(feature))
                    featureEnd = max(transProteinData.get(j).SeqEnd, featureEnd);
            }

            boolean pfPreserved = proteinFeaturePreserved(downTrans, true, featureStart, featureEnd);

            if(!pfPreserved && downTrans.gene().StableId.equals(fusion.upstreamTrans().gene().StableId))
            {
                // for same gene fusions, check whether the upstream transcript section preserves this feature
                pfPreserved = proteinFeaturePreserved(fusion.upstreamTrans(), false, featureStart, featureEnd);
            }

            downTrans.addProteinFeature(feature, pfPreserved);

            processedFeatures.add(feature);
        }

    }

    private static boolean proteinFeaturePreserved(
            final BreakendTransData transcript, boolean isDownstream, int featureStart, int featureEnd)
    {
        boolean featurePreserved;

        if(transcript.preCoding())
        {
            featurePreserved = isDownstream;
        }
        else if(transcript.postCoding())
        {
            featurePreserved = !isDownstream;
        }
        else
        {
            if (isDownstream)
            {
                // coding must start before the start of the feature for it to be preserved
                int featureCodingBaseStart = featureStart * 3;
                int svCodingBaseStart = transcript.CodingBases;
                featurePreserved = (featureCodingBaseStart >= svCodingBaseStart);
            }
            else
            {
                int svCodingBaseEnd = transcript.CodingBases;
                int featureCodingBaseEnd = featureEnd * 3;
                featurePreserved = (featureCodingBaseEnd <= svCodingBaseEnd);
            }
        }

        return featurePreserved;
    }

    private void setKnownFusionType(GeneFusion geneFusion)
    {
        final String upGene = geneFusion.transcripts()[FS_UP].gene().GeneName;
        final String downGene = geneFusion.transcripts()[FS_DOWN].gene().GeneName;

        if(mKnownFusionCache.hasKnownFusion(upGene, downGene))
        {
            geneFusion.setKnownType(KNOWN_PAIR);
            return;
        }

        if(upGene.equals(downGene))
        {
            if(mKnownFusionCache.withinKnownExonRanges(
                    EXON_DEL_DUP, geneFusion.transcripts()[FS_UP].transName(),
                    geneFusion.getBreakendExon(true), geneFusion.getFusedExon(true),
                    geneFusion.getBreakendExon(false), geneFusion.getFusedExon(false)))
            {
                geneFusion.setKnownType(EXON_DEL_DUP);
                geneFusion.setKnownExons();
                return;
            }

            // cannot be anything else, including promiscuous
            geneFusion.setKnownType(NONE);
            return;
        }

        if(mKnownFusionCache.hasPromiscuousThreeGene(downGene))
        {
            geneFusion.setKnownType(PROMISCUOUS_3);
            geneFusion.isPromiscuous()[FS_DOWN] = true;

            if(mKnownFusionCache.withinPromiscuousExonRange(
                    PROMISCUOUS_3, geneFusion.transcripts()[FS_DOWN].transName(),
                    geneFusion.getBreakendExon(false), geneFusion.getFusedExon(false)))
            {
                geneFusion.setKnownExons();
            }
        }

        if(mKnownFusionCache.hasPromiscuousFiveGene(upGene))
        {
            // will override promiscuous 3 but will show as PROM_BOTH in subsequent output
            geneFusion.setKnownType(PROMISCUOUS_5);
            geneFusion.isPromiscuous()[FS_UP] = true;

            if(mKnownFusionCache.withinPromiscuousExonRange(
                    PROMISCUOUS_5, geneFusion.transcripts()[FS_UP].transName(),
                    geneFusion.getBreakendExon(true), geneFusion.getFusedExon(true)))
            {
                geneFusion.setKnownExons();
            }
        }

        if(geneFusion.knownType() == PROMISCUOUS_5 || geneFusion.knownType() == PROMISCUOUS_3)
        {
            if(mKnownFusionCache.isHighImpactPromiscuous(geneFusion.knownType(), geneFusion.geneName(FS_UP), geneFusion.geneName(FS_DOWN)))
            {
                geneFusion.setHighImpactPromiscuous();
            }
        }
    }
}
