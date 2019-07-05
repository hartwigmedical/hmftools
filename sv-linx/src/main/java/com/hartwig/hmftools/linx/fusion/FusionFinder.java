package com.hartwig.hmftools.linx.fusion;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_3P_PROM;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_5P_PROM;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_BOTH_PROM;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_KNOWN;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_NONE;

import static com.hartwig.hmftools.linx.fusion.KnownFusionData.FUSION_PAIRS_CSV;
import static com.hartwig.hmftools.linx.fusion.KnownFusionData.PROMISCUOUS_FIVE_CSV;
import static com.hartwig.hmftools.linx.fusion.KnownFusionData.PROMISCUOUS_THREE_CSV;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.FusionAnnotations;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptProteinData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class FusionFinder
{

    private static final int EXON_THRESHOLD = 1;

    private KnownFusionData mKnownFusionData;
    boolean mHasValidConfigData;

    private SvGeneTranscriptCollection mGeneTranscriptCollection;
    private List<String> mProteinsRequiredKept;
    private List<String> mProteinsRequiredLost;

    private final String mOutputDir;
    private BufferedWriter mFusionWriter;
    private static boolean mLogInvalidReasons;

    private static final Logger LOGGER = LogManager.getLogger(FusionFinder.class);

    public FusionFinder(final CommandLine cmd, final SvGeneTranscriptCollection geneTransCache, final String outputDir)
    {
        mGeneTranscriptCollection = geneTransCache;
        mOutputDir = outputDir;

        mKnownFusionData = null;
        mHasValidConfigData = false;

        mFusionWriter = null;

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
        if(cmd.hasOption(FUSION_PAIRS_CSV) || cmd.hasOption(PROMISCUOUS_FIVE_CSV) || cmd.hasOption(PROMISCUOUS_THREE_CSV))
        {
            mKnownFusionData = new KnownFusionData();

            if(mKnownFusionData.loadFromFile(cmd))
            {
                LOGGER.debug("loaded known fusion data");
                mHasValidConfigData = true;
            }
        }
    }

    public void setHasValidConfigData(boolean toggle) { mHasValidConfigData = toggle; }
    public void setLogInvalidReasons(boolean toggle) { mLogInvalidReasons = toggle; }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(FUSION_PAIRS_CSV, true, "Path towards a CSV containing white-listed gene fusion pairs.");
        options.addOption(PROMISCUOUS_FIVE_CSV, true, "Path towards a CSV containing white-listed promiscuous 5' genes.");
        options.addOption(PROMISCUOUS_THREE_CSV, true, "Path towards a CSV containing white-listed promiscuous 3' genes.");
    }

    public final KnownFusionData getKnownFusionDatal() { return mKnownFusionData; }

    public static final String INVALID_REASON_ORIENTATION = "Orientation";
    public static final String INVALID_REASON_PHASING = "Unphased";
    public static final String INVALID_REASON_CODING_TYPE = "Coding";

    public final List<GeneFusion> findFusions(
            final List<GeneAnnotation> breakendGenes1, final List<GeneAnnotation> breakendGenes2,
            boolean requirePhaseMatch, boolean allowExonSkipping, @Nullable  List<String> invalidReasons, boolean setReportable)
    {
        final List<GeneFusion> potentialFusions = Lists.newArrayList();

        if(!mHasValidConfigData)
            return potentialFusions;

        for (final GeneAnnotation startGene : breakendGenes1)
        {
            // left is upstream, right is downstream
            boolean startUpstream = startGene.isUpstream();

            for (final GeneAnnotation endGene : breakendGenes2)
            {
                boolean endUpstream = endGene.isUpstream();

                if (startUpstream == endUpstream)
                {
                    if(invalidReasons!= null && !invalidReasons.contains(INVALID_REASON_ORIENTATION))
                        invalidReasons.add(INVALID_REASON_ORIENTATION);

                    continue;
                }

                for (final Transcript startTrans : startGene.transcripts())
                {
                    for (final Transcript endTrans : endGene.transcripts())
                    {
                        final Transcript upstreamTrans = startUpstream ? startTrans : endTrans;
                        final Transcript downstreamTrans = !startUpstream ? startTrans : endTrans;

                        GeneFusion geneFusion = checkFusionLogic(upstreamTrans, downstreamTrans, requirePhaseMatch,
                                allowExonSkipping, invalidReasons);

                        if(geneFusion == null)
                            continue;

                        geneFusion.setKnownFusionType(getKnownFusionType(upstreamTrans, downstreamTrans));

                        potentialFusions.add(geneFusion);
                    }
                }
            }
        }

        if(setReportable)
            setReportableGeneFusions(potentialFusions);

        return potentialFusions;
    }

    private static void logInvalidReasonInfo(final String reason, final Transcript trans1, final Transcript trans2,
            final String reasonType, @Nullable List<String> invalidReasons)
    {
        if(invalidReasons != null && !invalidReasons.contains(reasonType))
            invalidReasons.add(reasonType);

        if(!mLogInvalidReasons)
            return;

        if(trans2 == null)
            LOGGER.debug("transcript({}:{}) invalid({})", trans1.geneName(), trans1.StableId, reason);
        else
            LOGGER.debug("transcripts({}:{} and {}:{}) invalid({})",
                    trans1.geneName(), trans1.StableId, trans2.geneName(), trans2.StableId, reason);
    }

    public static boolean validFusionTranscript(final Transcript transcript)
    {
        // check any conditions which would preclude this transcript being a part of a fusion no matter the other end
        if(transcript.postCoding())
            return false;

        if(transcript.isUpstream())
        {
            if(transcript.isPromoter())
                return false;

            if(!transcript.isDisruptive())
                return false;
        }
        else
        {
            if(transcript.nonCoding())
                return false;

            if(transcript.ExonMax == 1)
                return false;
        }

        return true;
    }

    public static GeneFusion checkFusionLogic(final Transcript upstreamTrans, final Transcript downstreamTrans)
    {
        return checkFusionLogic(upstreamTrans, downstreamTrans, true, false, null);
    }

    public static GeneFusion checkFusionLogic(final Transcript upstreamTrans, final Transcript downstreamTrans,
            boolean requirePhaseMatch, boolean allowExonSkipping, @Nullable List<String> invalidReasons)
    {
        // see SV Fusions document for permitted combinations
        boolean checkExactMatch = false;

        if(!validFusionTranscript(upstreamTrans) || !validFusionTranscript(downstreamTrans))
        {
            logInvalidReasonInfo("invalid trans", upstreamTrans, downstreamTrans, INVALID_REASON_CODING_TYPE, invalidReasons);
            return null;
        }

        if(upstreamTrans.preCoding())
        {
            if(upstreamTrans.isExonic() && !downstreamTrans.isExonic())
            {
                logInvalidReasonInfo("precoding exonic to non-exonic", upstreamTrans, downstreamTrans,
                        INVALID_REASON_CODING_TYPE, invalidReasons);
                return null;
            }
            else if(downstreamTrans.isCoding())
            {
                logInvalidReasonInfo("pre-coding to coding", upstreamTrans, downstreamTrans,
                        INVALID_REASON_CODING_TYPE, invalidReasons);
                return null;
            }
        }
        else if(upstreamTrans.isCoding())
        {
            if(!downstreamTrans.isCoding())
            {
                logInvalidReasonInfo("coding to non-coding", upstreamTrans, downstreamTrans,
                        INVALID_REASON_CODING_TYPE, invalidReasons);
                return null;
            }

            if(upstreamTrans.isExonic())
            {
                if(!downstreamTrans.isExonic() && (!allowExonSkipping || !downstreamTrans.isIntronic()))
                {
                    logInvalidReasonInfo("coding exonic to non-exonic", upstreamTrans, downstreamTrans,
                            INVALID_REASON_CODING_TYPE, invalidReasons);
                    return null;
                }

                if(upstreamTrans.parent().id() != downstreamTrans.parent().id())
                {
                    logInvalidReasonInfo("up coding exonic diff SVs", upstreamTrans, downstreamTrans,
                            INVALID_REASON_CODING_TYPE, invalidReasons);
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
                logInvalidReasonInfo("up non-coding exonic to down non-exonic", upstreamTrans, downstreamTrans,
                        INVALID_REASON_CODING_TYPE, invalidReasons);
                return null;
            }
            else if(downstreamTrans.isCoding())
            {
                logInvalidReasonInfo("up non-coding to down-coding", upstreamTrans, downstreamTrans,
                        INVALID_REASON_CODING_TYPE, invalidReasons);
                return null;
            }
        }

        if (!isPotentiallyRelevantFusion(upstreamTrans, downstreamTrans))
        {
            logInvalidReasonInfo("irrelevant fusion", upstreamTrans, downstreamTrans, INVALID_REASON_CODING_TYPE, invalidReasons);
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
                logInvalidReasonInfo("downstream last exon", upstreamTrans, downstreamTrans,
                        INVALID_REASON_CODING_TYPE, invalidReasons);
                return null;
            }
        }

        if(checkExactMatch)
        {
            phaseMatched = exonToExonInPhase(upstreamTrans, downstreamTrans);

            if(!phaseMatched)
            {
                logInvalidReasonInfo("exact phase match", upstreamTrans, downstreamTrans, INVALID_REASON_PHASING, invalidReasons);
            }

            if(phaseMatched || !requirePhaseMatch)
            {
                return new GeneFusion(upstreamTrans, downstreamTrans, phaseMatched, true);
            }
        }
        else
        {
            // just check for a phasing match
            phaseMatched = upstreamTrans.ExonUpstreamPhase == downstreamTrans.ExonDownstreamPhase;

            if(!phaseMatched && allowExonSkipping)
            {
                // check for a match within the alterative phasings from upstream and downstream of the breakend
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

            if(!phaseMatched)
            {
                logInvalidReasonInfo("inexact unphased", upstreamTrans, downstreamTrans, INVALID_REASON_PHASING, invalidReasons);
            }

            if(phaseMatched || !requirePhaseMatch)
            {
                GeneFusion fusion = new GeneFusion(upstreamTrans, downstreamTrans, phaseMatched, true);
                fusion.setExonsSkipped(phaseExonsSkippedUp, phaseExonsSkippedDown);
                return fusion;
            }
        }

        return null;
    }

    private static boolean exonToExonInPhase(final Transcript transUp, final Transcript transDown)
    {
        // check phasing and offset since exon start or coding start
        transUp.setExonicCodingBase();
        transDown.setExonicCodingBase();

        return transUp.exactCodingBase() == transDown.exactCodingBase();
    }

    private static boolean isPotentiallyRelevantFusion(final Transcript t1, final Transcript t2)
    {
        if(!t1.geneName().equals(t2.geneName()))
            return true;

        // skip fusions between different transcripts in the same gene,
        if (!t1.StableId.equals(t2.StableId))
            return false;

        if(t1.nonCoding())
            return false;

        // skip fusions within the same intron
        if(t1.isIntronic() && t2.isIntronic() && t1.ExonUpstream == t2.ExonUpstream)
            return false;

        return true;
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

        if(intragenic(reportableFusion.upstreamTrans(), reportableFusion.downstreamTrans())
        && reportableFusion.upstreamTrans().ExonUpstreamPhase == -1)
        {
            return;
        }

        reportableFusion.setReportable(true);

        // check impact on protein regions
        setFusionProteinFeatures(reportableFusion);

        if (reportableFusion.getKnownFusionType() != REPORTABLE_TYPE_KNOWN)
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

        final List<TranscriptProteinData> transProteinData = mGeneTranscriptCollection.getTranscriptProteinDataMap().get(downTrans.TransId);

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

            addProteinFeature(downTrans, true, feature, featureStart, featureEnd);
            processedFeatures.add(feature);
        }

    }

    private void addProteinFeature(final Transcript transcript, boolean isDownstream, final String feature, int featureStart, int featureEnd)
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
                int codingBaseStart = transcript.totalCodingBases() - transcript.calcCodingBases(!isDownstream);
                int featureCodingBaseStart = featureStart * 3;
                featurePreserved = (codingBaseStart <= featureCodingBaseStart);
            }
            else
            {
                int codingBaseEnd = transcript.calcCodingBases(!isDownstream);
                int featureCodingBaseEnd = featureEnd * 3;
                featurePreserved = (featureCodingBaseEnd <= codingBaseEnd);
            }
        }

        transcript.addProteinFeature(feature, featurePreserved);
    }

    public static String TRANSCRIPT_PROTEIN_CODING = "protein_coding";
    public static String TRANSCRIPT_NONSENSE_MED_DECAY = "nonsense_mediated_decay";

    private static int MAX_UPSTREAM_DISTANCE_KNOWN = 100000;
    private static int MAX_UPSTREAM_DISTANCE_OTHER = 10000;

    public static boolean couldBeReportable(GeneFusion fusion)
    {
        if(!fusion.phaseMatched())
            return false;

        // first check whether a fusion is known or not - a key requirement of it being potentially reportable
        if (fusion.getKnownFusionType() == REPORTABLE_TYPE_NONE)
            return false;

        // set limits on how far upstream the breakend can be - adjusted for whether the fusions is known or not
        int maxUpstreamDistance = fusion.getKnownFusionType() == REPORTABLE_TYPE_KNOWN ?
                MAX_UPSTREAM_DISTANCE_KNOWN : MAX_UPSTREAM_DISTANCE_OTHER;

        final Transcript upTrans = fusion.upstreamTrans();
        final Transcript downTrans = fusion.downstreamTrans();

        if(upTrans.getDistanceUpstream() > maxUpstreamDistance || downTrans.getDistanceUpstream() > maxUpstreamDistance)
            return false;

        if(downTrans.bioType().equals(TRANSCRIPT_NONSENSE_MED_DECAY))
            return false;

        if(downTrans.exonDistanceUp() < 0)
            return false;

        if(fusion.getKnownFusionType() != REPORTABLE_TYPE_KNOWN
        && (fusion.getExonsSkipped(true) > 0 || fusion.getExonsSkipped(false) > 0))
            return false;

        return true;
    }

    public static GeneFusion determineReportableFusion(final List<GeneFusion> fusions, boolean requireReportable)
    {
        GeneFusion reportableFusion = null;

        // form a score by allocating 0/1 or length value to each power 10 descending
        long highestScore = 0;

        for(final GeneFusion fusion : fusions)
        {
            if(requireReportable && !couldBeReportable(fusion))
                continue;

            // first check whether a fusion is known or not - a key requirement of it being potentially reportable
            final Transcript upTrans = fusion.upstreamTrans();
            final Transcript downTrans = fusion.downstreamTrans();

            /* prioritisation rules:
            - prioritise fusions which don't skip exons
            - take both canonical if possible
            - favour 3' over 5' by canoncial, protein-coding then coding bases (or exon count if not coding)
            */

            long transScore = 0;
            long factor = 1000000000;

            if(fusion.getExonsSkipped(true) == 0 && fusion.getExonsSkipped(false) == 0)
                transScore += factor;

            factor /= 10;

            if(downTrans.isCanonical() && upTrans.isCanonical())
                transScore += factor;

            factor /= 10;

            if(downTrans.isCanonical())
                transScore += factor;

            factor /= 10;

            if(downTrans.bioType().equals(TRANSCRIPT_PROTEIN_CODING))
                transScore += factor;

            factor /= 100;

            long length = downTrans.isCoding() ? downTrans.calcCodingBases(false) : downTrans.ExonMax;

            // will be a range between 1-99 * current factor
            length = min(round(length/10), 99);
            transScore += length * factor;

            factor /= 10;

            if(upTrans.isCanonical())
                transScore += factor;

            factor /= 10;

            if(upTrans.bioType().equals(TRANSCRIPT_PROTEIN_CODING))
                transScore += factor;

            factor /= 100;

            length = upTrans.isCoding() ? upTrans.calcCodingBases(true) : upTrans.ExonMax;
            length = min(round(length/10), 99);
            transScore += length * factor;

            if(transScore > highestScore)
            {
                reportableFusion = fusion;
                highestScore = transScore;
            }
        }

        return reportableFusion;
    }

    public String getKnownFusionType(final Transcript upTrans, final Transcript downTrans)
    {
        if(mHasValidConfigData && mKnownFusionData == null)
            return REPORTABLE_TYPE_KNOWN;

        final String upGene = upTrans.parent().GeneName;
        final String downGene = downTrans.parent().GeneName;

        if(mKnownFusionData.hasKnownFusion(upGene, downGene))
            return REPORTABLE_TYPE_KNOWN;

        boolean fivePrimeMatch = mKnownFusionData.hasPromiscuousFiveGene(upGene);
        boolean threePrimeMatch = mKnownFusionData.hasPromiscuousThreeGene(downGene);

        boolean intergenicPromiscuousMatch = mKnownFusionData.intergenicPromiscuousMatch(upGene, downGene);

        boolean intragenicPromiscuousMatch = mKnownFusionData.intragenicPromiscuousMatch(upGene, downGene)
                && downTrans.ExonDownstream - upTrans.ExonUpstream > EXON_THRESHOLD;

        if(intergenicPromiscuousMatch || intragenicPromiscuousMatch)
        {
            if (fivePrimeMatch && threePrimeMatch)
                return REPORTABLE_TYPE_BOTH_PROM;
            else if (fivePrimeMatch)
                return REPORTABLE_TYPE_5P_PROM;
            else if (threePrimeMatch)
                return REPORTABLE_TYPE_3P_PROM;
            else
                return REPORTABLE_TYPE_NONE;
        }

        return REPORTABLE_TYPE_NONE;
    }

    private static boolean intragenic(final Transcript upstream, final Transcript downstream)
    {
        return upstream.parent().synonyms().stream().anyMatch(downstream.parent().synonyms()::contains);
    }

    public void initialiseOutputFile(final String fileName, String annotatationHeaders)
    {
        try
        {
            if(mFusionWriter == null)
            {
                String outputFilename = mOutputDir;

                if (!outputFilename.endsWith(File.separator))
                    outputFilename += File.separator;

                outputFilename += fileName;

                mFusionWriter = createBufferedWriter(outputFilename, false);

                mFusionWriter.write("SampleId,Reportable,KnownType");

                if(annotatationHeaders.isEmpty())
                    mFusionWriter.write(",ClusterId,ClusterCount,ClusterInfo");
                else
                    mFusionWriter.write(String.format(",%s", annotatationHeaders));


                mFusionWriter.write(",SvIdUp,ChrUp,PosUp,OrientUp,TypeUp,PloidyUp,GeneIdUp,GeneNameUp,ChrBandUp,TranscriptUp,StrandUp,RegionTypeUp,CodingTypeUp");
                mFusionWriter.write(",ExonUp,PhaseUp,ExonMaxUp,DisruptiveUp,ExactBaseUp,CodingBasesUp,TotalCodingUp");
                mFusionWriter.write(",CodingStartUp,CodingEndUp,TransStartUp,TransEndUp,DistancePrevUp,CanonicalUp,BiotypeUp");

                mFusionWriter.write(",SvIdDown,ChrDown,PosDown,OrientDown,TypeDown,PloidyDown,GeneIdDown,GeneNameDown,ChrBandDown,TranscriptDown,StrandDown,RegionTypeDown,CodingTypeDown");
                mFusionWriter.write(",ExonDown,PhaseDown,ExonMaxDown,DisruptiveDown,ExactBaseDown,CodingBasesDown,TotalCodingDown");
                mFusionWriter.write(",CodingStartDown,CodingEndDown,TransStartDown,TransEndDown,DistancePrevDown,CanonicalDown,BiotypeDown");

                mFusionWriter.write(",ProteinsKept,ProteinsLost,ExonsSkippedUp,ExonsSkippedDown");
                mFusionWriter.newLine();
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing fusions: {}", e.toString());
        }
    }

    public void writeFusionData(final GeneFusion fusion, final String sampleId)
    {
        if(mFusionWriter == null)
            return;

        try
        {
            BufferedWriter writer = mFusionWriter;

            final Transcript upTrans = fusion.upstreamTrans();
            final Transcript downTrans = fusion.downstreamTrans();

            final GeneAnnotation upGene = upTrans.parent();
            final GeneAnnotation downGene = downTrans.parent();

            final FusionAnnotations annotations = fusion.getAnnotations();
            String annotationsStr = ",,";

            if(annotations != null)
            {
                // PhaseMatched,ClusterId,ClusterCount,ResolvedType,OverlapUp,OverlapDown,ChainInfo
                annotationsStr = String.format("%s,%d,%d,%s",
                        fusion.phaseMatched(), annotations.clusterId(), annotations.clusterCount(), annotations.resolvedType());

                String defaultValues = ",0:false;0;0;0;false";
                if(annotations.disruptionUp() != null)
                {
                    annotationsStr += String.format(",%d;%s;%d;%d;%d;%s",
                            annotations.disruptionUp().facingBreakends(), annotations.disruptionUp().allLinksAssembled(),
                            annotations.disruptionUp().totalBreakends(), annotations.disruptionUp().minDistance(),
                            annotations.disruptionUp().disruptedExons(), annotations.disruptionUp().transcriptTerminated());
                }
                else
                {
                    annotationsStr += defaultValues;

                }

                if(annotations.disruptionDown() != null)
                {
                    annotationsStr += String.format(",%d;%s;%d;%d;%d;%s",
                            annotations.disruptionDown().facingBreakends(), annotations.disruptionDown().allLinksAssembled(),
                            annotations.disruptionDown().totalBreakends(), annotations.disruptionDown().minDistance(),
                            annotations.disruptionDown().disruptedExons(), annotations.disruptionDown().transcriptTerminated());
                }
                else
                {
                    annotationsStr += defaultValues;
                }

                if(annotations.chainInfo() != null)
                {
                    annotationsStr += String.format(",%d;%d;%d;%s;%s",
                            annotations.chainInfo().chainId(), annotations.chainInfo().chainLinks(), annotations.chainInfo().chainLength(),
                            annotations.chainInfo().validTraversal(), annotations.chainInfo().traversalAssembled());
                }
                else
                {
                    annotationsStr += ",-1;0;0;true;false";
                }
            }

            writer.write(String.format("%s,%s,%s,%s",
                    sampleId, fusion.reportable(), fusion.getKnownFusionType(), annotationsStr));

            // write upstream SV, transcript and exon info
            writer.write(
                    String.format(",%d,%s,%d,%d,%s,%.6f",
                            upGene.id(), upGene.chromosome(), upGene.position(), upGene.orientation(),
                            upGene.type(), upGene.ploidy()));

            writer.write(
                    String.format(",%s,%s,%s,%s,%d,%s,%s",
                            upGene.StableId, upGene.GeneName, upGene.karyotypeBand(), upTrans.StableId,
                            upGene.Strand, upTrans.regionType(), upTrans.codingType()));

            writer.write(
                    String.format(",%d,%d,%d,%s",
                            upTrans.ExonUpstream, upTrans.ExonUpstreamPhase, upTrans.ExonMax, upTrans.isDisruptive()));
            writer.write(
                    String.format(",%d,%d,%d,%d,%d,%d,%d,%d,%s,%s",
                            upTrans.exactCodingBase(), upTrans.calcCodingBases(true), upTrans.totalCodingBases(),
                            upTrans.codingStart(), upTrans.codingEnd(), upTrans.TranscriptStart, upTrans.TranscriptEnd,
                            upTrans.exonDistanceUp(), upTrans.isCanonical(), upTrans.bioType()));

            writer.write(
                    String.format(",%d,%s,%d,%d,%s,%.6f",
                            downGene.id(), downGene.chromosome(), downGene.position(), downGene.orientation(),
                            downGene.type(), downGene.ploidy()));

            writer.write(
                    String.format(",%s,%s,%s,%s,%d,%s,%s",
                            downGene.StableId, downGene.GeneName, downGene.karyotypeBand(), downTrans.StableId,
                            downGene.Strand, downTrans.regionType(), downTrans.codingType()));

            writer.write(
                    String.format(",%d,%d,%d,%s",
                            downTrans.ExonDownstream, downTrans.ExonDownstreamPhase, downTrans.ExonMax, downTrans.isDisruptive()));

            writer.write(
                    String.format(",%d,%d,%d,%d,%d,%d,%d,%d,%s,%s,%s,%s",
                            downTrans.exactCodingBase(), downTrans.calcCodingBases(false), downTrans.totalCodingBases(),
                            downTrans.codingStart(), downTrans.codingEnd(), downTrans.TranscriptStart, downTrans.TranscriptEnd,
                            downTrans.exonDistanceUp(), downTrans.isCanonical(), downTrans.bioType(),
                            downTrans.getProteinFeaturesKept(), downTrans.getProteinFeaturesLost()));

            // move to phasing section once move to SVA fusions
            writer.write(String.format(",%d,%d", fusion.getExonsSkipped(true), fusion.getExonsSkipped(false)));

            writer.newLine();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing fusions: {}", e.toString());
        }
    }

    private static int MAX_PROMOTOR_DISTANCE_UP = 100000;

    public boolean isTranscriptBreakendViableForRnaBoundary(final Transcript trans, boolean isUpstream, long breakendPosition,
            long rnaPosition, boolean exactRnaPosition)
    {
        final TranscriptData transData = mGeneTranscriptCollection.getTranscriptData(trans.parent().StableId, trans.StableId);

        if (transData == null || transData.exons().isEmpty())
            return false;

        int strand = trans.parent().Strand;

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
                    // if the RNA boundary is at or before the 2nd exon (which has the first splice acceptor), then the breakend can
                    // be upstream as far the previous gene or 100K
                    int distanceUp = trans.exonDistanceUp();
                    long breakendDistance = abs(breakendPosition - rnaPosition);

                    if(breakendDistance > MAX_PROMOTOR_DISTANCE_UP || distanceUp < 0)
                        return false;
                    else
                        return true;
                }

                // breakend must fall at or before the RNA boundary but not further upstream than the previous splice acceptor
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

    public void onCompleted()
    {
        closeBufferedWriter(mFusionWriter);
    }
}
