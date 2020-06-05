package com.hartwig.hmftools.linx.fusion;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.ensemblcache.TranscriptProteinData.BIOTYPE_NONSENSE_MED_DECAY;
import static com.hartwig.hmftools.common.ensemblcache.TranscriptProteinData.BIOTYPE_PROTEIN_CODING;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWNSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UPSTREAM;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.EXON_DEL_DUP;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_PROMISCUOUS;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.NONE;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_3;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_5;

import static com.hartwig.hmftools.common.fusion.KnownFusionCache.FUSION_PAIRS_CSV;
import static com.hartwig.hmftools.common.fusion.KnownFusionCache.KNOWN_FUSIONS_FILE;
import static com.hartwig.hmftools.common.fusion.KnownFusionCache.PROMISCUOUS_FIVE_CSV;
import static com.hartwig.hmftools.common.fusion.KnownFusionCache.PROMISCUOUS_THREE_CSV;
import static com.hartwig.hmftools.common.fusion.Transcript.TRANS_REGION_TYPE_UPSTREAM;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.MAX_UPSTREAM_DISTANCE_KNOWN;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.MAX_UPSTREAM_DISTANCE_OTHER;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.REQUIRED_BIOTYPES;
import static com.hartwig.hmftools.linx.fusion.FusionReportability.checkProteinDomains;
import static com.hartwig.hmftools.linx.fusion.FusionReportability.couldBeReportable;
import static com.hartwig.hmftools.linx.fusion.FusionReportability.determineReportableFusion;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.fusion.GeneAnnotation;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.fusion.Transcript;
import com.hartwig.hmftools.common.ensemblcache.TranscriptProteinData;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.jetbrains.annotations.NotNull;

public class FusionFinder
{
    private final KnownFusionCache mKnownFusionCache;
    private boolean mHasValidConfigData;

    private EnsemblDataCache mGeneTransCache;
    private List<String> mProteinsRequiredKept;
    private List<String> mProteinsRequiredLost;

    private static boolean mLogInvalidReasons;

    private static final int EXON_THRESHOLD = 1;

    public FusionFinder(final CommandLine cmd, final EnsemblDataCache geneTransCache)
    {
        mGeneTransCache = geneTransCache;

        mKnownFusionCache = new KnownFusionCache();
        mHasValidConfigData = false;

        mProteinsRequiredKept = Lists.newArrayList();
        mProteinsRequiredLost = Lists.newArrayList();
        setRequiredProteins();

        if(cmd != null)
        {
            initialise(cmd);
        }

        mLogInvalidReasons = false;
    }

    private void initialise(@NotNull final CommandLine cmd)
    {
        if(mKnownFusionCache.loadFromFile(cmd))
        {
            LNX_LOGGER.debug("loaded known fusion data");
            mHasValidConfigData = true;
        }
    }

    public void setHasValidConfigData(boolean toggle) { mHasValidConfigData = toggle; }
    public boolean hasValidConfigData() { return mHasValidConfigData; }
    public void setLogInvalidReasons(boolean toggle) { mLogInvalidReasons = toggle; }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(KNOWN_FUSIONS_FILE, true, "Known fusion file");
        options.addOption(FUSION_PAIRS_CSV, true, "Known gene fusion pairs");
        options.addOption(PROMISCUOUS_FIVE_CSV, true, "Promiscuous 5' genes");
        options.addOption(PROMISCUOUS_THREE_CSV, true, "Promiscuous 3' genes");
    }

    public final KnownFusionCache getKnownFusionCache() { return mKnownFusionCache; }

    public static final String INVALID_REASON_ORIENTATION = "Orientation";
    public static final String INVALID_REASON_PHASING = "Unphased";
    public static final String INVALID_REASON_CODING_TYPE = "Coding";

    public final List<GeneFusion> findFusions(
            final List<GeneAnnotation> breakendGenes1, final List<GeneAnnotation> breakendGenes2,
            final FusionParameters params, boolean setReportable)
    {
        final List<GeneFusion> potentialFusions = Lists.newArrayList();

        if(!mHasValidConfigData)
            return potentialFusions;

        for (final GeneAnnotation startGene : breakendGenes1)
        {
            // left is upstream, right is downstream
            boolean startUpstream = startGene.isUpstream();
            boolean startIsIgRegion = mKnownFusionCache.withinIgRegion(startGene.chromosome(), startGene.position());

            for (final GeneAnnotation endGene : breakendGenes2)
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

                final GeneAnnotation upGene = startUpstream ? startGene : endGene;
                final GeneAnnotation downGene = !startUpstream ? startGene : endGene;

                boolean knownPair = mKnownFusionCache.hasKnownFusion(upGene.GeneName, downGene.GeneName);

                for(final Transcript upstreamTrans : upGene.transcripts())
                {
                    for(final Transcript downstreamTrans : downGene.transcripts())
                    {
                        GeneFusion geneFusion = checkFusionLogic(upstreamTrans, downstreamTrans, params, !knownPair);

                        if(geneFusion == null)
                            continue;

                        setKnownFusionType(geneFusion);
                        potentialFusions.add(geneFusion);
                    }
                }
            }
        }

        if(setReportable)
            setReportableGeneFusions(potentialFusions);

        return potentialFusions;
    }

    private void checkIgFusion(final GeneAnnotation startGene, final GeneAnnotation endGene, final List<GeneFusion> potentialFusions)
    {
        /* Criteria:
            - These are allowed to fuse with 5’UTR splice donors from 10kb upstream. Phasing is assumed to be -1.
            - In the case of IGH-BCL2 specifically, allow fusions in the 3’UTR region and up to 40k bases downstream of BCL2 ((common in Folicular Lymphomas[PP1] )
        */

        boolean startIsIgGene = mKnownFusionCache.matchesIgGene(startGene.chromosome(), startGene.position(), startGene.orientation());
        boolean endIsIgGene = !startIsIgGene && mKnownFusionCache.matchesIgGene(endGene.chromosome(), endGene.position(), endGene.orientation());

        if(!startIsIgGene && !endIsIgGene)
            return;

        final GeneAnnotation igGene = startIsIgGene ? startGene : endGene;
        final GeneAnnotation downGene = startIsIgGene ? endGene : startGene;

        KnownFusionType knownType = NONE;

        final List<Transcript> candidateTranscripts = downGene.transcripts().stream()
                .filter(x -> x.regionType().equals(TRANS_REGION_TYPE_UPSTREAM)).collect(Collectors.toList());

        KnownFusionData knownFusionData = mKnownFusionCache.getDataByType(IG_KNOWN_PAIR).stream()
                .filter(x -> x.ThreeGene.equals(downGene.GeneName))
                .filter(x -> x.withinIgRegion(igGene.chromosome(), igGene.position()))
                .findFirst().orElse(null);

        if(knownFusionData != null)
        {
            // a known IG-partner gene
            if(knownFusionData.igDownstreamDistance() > 0)
            {
                candidateTranscripts.addAll(downGene.transcripts().stream().filter(x -> x.postCoding()).collect(Collectors.toList()));
            }

            knownType = IG_KNOWN_PAIR;
        }
        else
        {
            knownFusionData = mKnownFusionCache.getDataByType(IG_PROMISCUOUS).stream()
                    .filter(x -> x.withinIgRegion(igGene.chromosome(), igGene.position()))
                    .findFirst().orElse(null);

            // check within the promiscuous region bounds
            if(knownFusionData == null)
                return;

            knownType = IG_PROMISCUOUS;
        }

        final Transcript upTrans = generateIgTranscript(igGene, knownFusionData);

        if(!candidateTranscripts.isEmpty())
        {
            for(final Transcript downTrans : candidateTranscripts)
            {
                GeneFusion fusion = new GeneFusion(upTrans, downTrans, true);
                fusion.setKnownType(knownType);
                potentialFusions.add(fusion);
            }
        }
    }

    private Transcript generateIgTranscript(final GeneAnnotation gene, final KnownFusionData knownFusionData)
    {
        return new Transcript(
                gene, 0, knownFusionData.FiveGene, 1,-1, 1, -1,
        0, 0,   1, true, knownFusionData.igRegion()[SE_START], knownFusionData.igRegion()[SE_END],
        null, null);
    }

    private static void logInvalidReasonInfo(final Transcript trans1, final Transcript trans2, final String reasonType, final String reason)
    {
        if(!mLogInvalidReasons)
            return;

        if(trans2 == null)
            LNX_LOGGER.trace("transcript({}:{}) invalid({}: {})", trans1.geneName(), trans1.StableId, reasonType, reason);
        else
            LNX_LOGGER.trace("transcripts({}:{} and {}:{}) invalid({}: {})",
                    trans1.geneName(), trans1.StableId, trans2.geneName(), trans2.StableId, reasonType, reason);
    }

    public static boolean validFusionTranscript(final Transcript transcript)
    {
        return validFusionTranscript(transcript, true, false);
    }

    private static boolean validFusionTranscript(
            final Transcript transcript, boolean requireUpstreamDisruptive, boolean requireUpstreamBiotypes)
    {
        // check any conditions which would preclude this transcript being a part of a fusion no matter the other end
        if(transcript.isUpstream())
        {
            if(transcript.isPromoter())
                return false;

            if(requireUpstreamDisruptive && !transcript.isDisruptive())
                return false;

            if(requireUpstreamBiotypes && !REQUIRED_BIOTYPES.contains(transcript.bioType()))
                return false;
        }
        else
        {
            if(transcript.postCoding())
                return false;

            if(transcript.nonCoding())
                return false;

            if(transcript.ExonMax == 1)
                return false;
        }

        return true;
    }

    public static GeneFusion checkFusionLogic(final Transcript upstreamTrans, final Transcript downstreamTrans, final FusionParameters params)
    {
        return checkFusionLogic(upstreamTrans, downstreamTrans, params, true);
    }

    private static GeneFusion checkFusionLogic(
            final Transcript upstreamTrans, final Transcript downstreamTrans, final FusionParameters params, boolean requireUpstreamDisruptive)
    {
        // see SV Fusions document for permitted combinations
        boolean checkExactMatch = false;

        if(!validFusionTranscript(upstreamTrans, requireUpstreamDisruptive, params.RequireUpstreamBiotypes)
        || !validFusionTranscript(downstreamTrans))
        {
            logInvalidReasonInfo(upstreamTrans, downstreamTrans, INVALID_REASON_CODING_TYPE, "invalid trans");
            return null;
        }

        if(upstreamTrans.preCoding())
        {
            if(upstreamTrans.isExonic() && !downstreamTrans.isExonic())
            {
                logInvalidReasonInfo(upstreamTrans, downstreamTrans, INVALID_REASON_CODING_TYPE, "precoding exonic to non-exonic");
                return null;
            }
            else if(downstreamTrans.isCoding())
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
            if(!downstreamTrans.isCoding())
            {
                logInvalidReasonInfo(upstreamTrans, downstreamTrans, INVALID_REASON_CODING_TYPE, "coding to non-coding");
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
            if(upstreamTrans.isExonic() && !downstreamTrans.isExonic())
            {
                logInvalidReasonInfo(upstreamTrans, downstreamTrans, INVALID_REASON_CODING_TYPE, "up non-coding exonic to down non-exonic");
                return null;
            }
            else if(downstreamTrans.isCoding())
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
            if(downstreamTrans.isExonic() && downstreamTrans.ExonDownstream == downstreamTrans.ExonMax && !downstreamTrans.preCoding())
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
                phaseMatched = upstreamTrans.ExonUpstreamPhase == downstreamTrans.ExonDownstreamPhase;
            }

            if(!phaseMatched && params.AllowExonSkipping && !upstreamTrans.gene().StableId.equals(downstreamTrans.gene().StableId))
            {
                // check for a match within the alternative phasings from upstream and downstream of the breakend
                for (Map.Entry<Integer, Integer> altPhasing : upstreamTrans.getAlternativePhasing().entrySet())
                {
                    if (altPhasing.getKey() == downstreamTrans.ExonDownstreamPhase)
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
                        if (altPhasing.getKey() == upstreamTrans.ExonUpstreamPhase)
                        {
                            phaseMatched = true;
                            phaseExonsSkippedDown = altPhasing.getValue();
                            break;
                        }
                    }
                }
            }

            if(upstreamTrans.isExonic() && !downstreamTrans.isExonic())
            {
                if(phaseExonsSkippedUp == 0)
                    phaseExonsSkippedUp = 1;
            }
            else if(!upstreamTrans.isExonic() && downstreamTrans.isExonic())
            {
                if(phaseExonsSkippedDown == 0)
                    phaseExonsSkippedDown = 1;
            }

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

    private static boolean exonToExonInPhase(final Transcript upTrans, final Transcript downTrans)
    {
        // check phasing and offset since exon start or coding start
        upTrans.setExonicCodingBase();
        downTrans.setExonicCodingBase();

        int upPhase = upTrans.exonicBasePhase();
        int downPhase = downTrans.exonicBasePhase();

        if(upPhase == -1 && downPhase == -1)
            return true;

        return ((upPhase + 1) % 3) == (downPhase % 3);
    }

    public static boolean isIrrelevantSameGene(final Transcript upTrans, final Transcript downTrans)
    {
        if(!upTrans.geneName().equals(downTrans.geneName()))
            return false;

        // skip fusions between different transcripts in the same gene,
        if (!upTrans.StableId.equals(downTrans.StableId))
            return true;

        if(upTrans.nonCoding())
            return true;

        // skip fusions within the same intron
        if(upTrans.isIntronic() && downTrans.isIntronic() && upTrans.ExonUpstream == downTrans.ExonUpstream)
            return true;

        return false;
    }

    private void setRequiredProteins()
    {
        mProteinsRequiredLost.add("Raf-like Ras-binding");

        mProteinsRequiredKept.add("Ets domain");
        mProteinsRequiredKept.add("Protein kinase domain");
        mProteinsRequiredKept.add("Epidermal growth factor-like domain");
        mProteinsRequiredKept.add("Ankyrin repeat-containing domain");
        mProteinsRequiredKept.add("Basic-leucine zipper domain");
        mProteinsRequiredKept.add("High mobility group box domain");
        mProteinsRequiredKept.add("Bromodomain");
    }

    public void setReportableGeneFusions(final List<GeneFusion> fusions)
    {
        if(fusions.isEmpty())
            return;

        GeneFusion reportableFusion = determineReportableFusion(fusions, true);

        if(reportableFusion == null)
            return;

        reportableFusion.setReportable(true);

        // check impact on protein regions
        setFusionProteinFeatures(reportableFusion);

        if (checkProteinDomains(reportableFusion.knownType()))
        {
            final Transcript downTrans = reportableFusion.downstreamTrans();
            long requiredKeptButLost = mProteinsRequiredKept.stream().filter(f -> downTrans.getProteinFeaturesLost().contains(f)).count();
            long requiredLostButKept = mProteinsRequiredLost.stream().filter(f -> downTrans.getProteinFeaturesKept().contains(f)).count();

            if (requiredKeptButLost > 0 || requiredLostButKept > 0)
                reportableFusion.setReportable(false);
        }
    }

    public void setFusionProteinFeatures(GeneFusion fusion)
    {
        final Transcript downTrans = fusion.downstreamTrans();

        if(!downTrans.getProteinFeaturesKept().isEmpty() || !downTrans.getProteinFeaturesLost().isEmpty())
            return;

        if(downTrans.nonCoding())
            return;

        final List<TranscriptProteinData> transProteinData = mGeneTransCache.getTranscriptProteinDataMap().get(downTrans.TransId);

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

            boolean pfPreserved = proteinFeaturePreserved(downTrans, true, feature, featureStart, featureEnd);

            if(!pfPreserved && downTrans.gene().StableId.equals(fusion.upstreamTrans().gene().StableId))
            {
                // for same gene fusions, check whether the upstream transcript section preserves this feature
                pfPreserved = proteinFeaturePreserved(fusion.upstreamTrans(), false, feature, featureStart, featureEnd);
            }

            downTrans.addProteinFeature(feature, pfPreserved);

            processedFeatures.add(feature);
        }

    }

    private static boolean proteinFeaturePreserved(
            final Transcript transcript, boolean isDownstream, final String feature, int featureStart, int featureEnd)
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
                int svCodingBaseStart = transcript.codingBases();
                featurePreserved = (featureCodingBaseStart >= svCodingBaseStart);
            }
            else
            {
                int svCodingBaseEnd = transcript.codingBases();
                int featureCodingBaseEnd = featureEnd * 3;
                featurePreserved = (featureCodingBaseEnd <= svCodingBaseEnd);
            }
        }

        return featurePreserved;
    }

    private void setKnownFusionType(GeneFusion geneFusion)
    {
        final String upGene = geneFusion.transcripts()[FS_UPSTREAM].gene().GeneName;
        final String downGene = geneFusion.transcripts()[FS_DOWNSTREAM].gene().GeneName;

        if(mKnownFusionCache.hasKnownFusion(upGene, downGene))
        {
            geneFusion.setKnownType(KNOWN_PAIR);
            return;
        }

        if(upGene.equals(downGene))
        {
            if(mKnownFusionCache.isExonDelDup(
                    upGene, geneFusion.transcripts()[FS_UPSTREAM].StableId,
                    geneFusion.getFusedExon(true), geneFusion.getFusedExon(false)))
            {
                geneFusion.setKnownType(EXON_DEL_DUP);
                return;
            }

            // cannot be anything else, including promiscuous
            geneFusion.setKnownType(NONE);
            return;
        }

        if(mKnownFusionCache.hasPromiscuousThreeGene(downGene))
        {
            geneFusion.setKnownType(PROMISCUOUS_3);
            geneFusion.isPromiscuous()[FS_DOWNSTREAM] = true;
        }

        if(mKnownFusionCache.hasPromiscuousFiveGene(upGene))
        {
            // will override promiscuous 3 but will show as PROM_BOTH in subsequent output
            geneFusion.setKnownType(PROMISCUOUS_5);
            geneFusion.isPromiscuous()[FS_UPSTREAM] = true;
        }
    }
}
