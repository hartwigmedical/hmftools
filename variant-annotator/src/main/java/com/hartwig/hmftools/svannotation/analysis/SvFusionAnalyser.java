package com.hartwig.hmftools.svannotation.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation.isUpstream;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_3P_PROM;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_5P_PROM;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_BOTH_PROM;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_KNOWN;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_NONE;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.EXON_RANK_MAX;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.EXON_RANK_MIN;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptProteinData;
import com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvFusionAnalyser
{

    public static final String FUSION_PAIRS_CSV = "fusion_pairs_csv";
    public static final String PROMISCUOUS_FIVE_CSV = "promiscuous_five_csv";
    public static final String PROMISCUOUS_THREE_CSV = "promiscuous_three_csv";

    private static final int EXON_THRESHOLD = 1;

    private KnownFusionsModel mKnownFusionsModel;

    private SvGeneTranscriptCollection mGeneTranscriptCollection;
    private List<String> mProteinsRequiredKept;
    private List<String> mProteinsRequiredLost;

    private final String mOutputDir;
    private BufferedWriter mFusionWriter;

    private static final Logger LOGGER = LogManager.getLogger(SvFusionAnalyser.class);

    public SvFusionAnalyser(final CommandLine cmd, final SvGeneTranscriptCollection geneTranscriptCollection, final String outputDir)
    {
        mGeneTranscriptCollection = geneTranscriptCollection;
        mOutputDir = outputDir;

        mKnownFusionsModel = null;

        mFusionWriter = null;

        mProteinsRequiredKept = Lists.newArrayList();
        mProteinsRequiredLost = Lists.newArrayList();
        setRequiredProteins();

        if(cmd != null)
        {
            loadConfig(cmd);
        }
    }

    private void loadConfig(final CommandLine cmd)
    {
        if(cmd.hasOption(FUSION_PAIRS_CSV) && cmd.hasOption(PROMISCUOUS_FIVE_CSV) && cmd.hasOption(PROMISCUOUS_THREE_CSV))
        {
            try
            {
                mKnownFusionsModel = KnownFusionsModel.fromInputStreams(new FileInputStream(cmd.getOptionValue(FUSION_PAIRS_CSV)),
                        new FileInputStream(cmd.getOptionValue(PROMISCUOUS_FIVE_CSV)),
                        new FileInputStream(cmd.getOptionValue(PROMISCUOUS_THREE_CSV)));

                LOGGER.debug("loaded known fusion data");
            }
            catch (IOException e)
            {
                LOGGER.warn("no known fusion files loaded");
            }
        }
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(FUSION_PAIRS_CSV, true, "Path towards a CSV containing white-listed gene fusion pairs.");
        options.addOption(PROMISCUOUS_FIVE_CSV, true, "Path towards a CSV containing white-listed promiscuous 5' genes.");
        options.addOption(PROMISCUOUS_THREE_CSV, true, "Path towards a CSV containing white-listed promiscuous 3' genes.");
    }

    public final List<GeneFusion> findFusions(final List<StructuralVariantAnnotation> annotations)
    {
        LOGGER.debug("finding fusions in {} annotations", annotations.size());

        List<GeneFusion> fusions = Lists.newArrayList();

        for (final StructuralVariantAnnotation annotation : annotations)
        {
            List<GeneFusion> svFusions = findFusions(annotation.start(), annotation.end());

            fusions.addAll(svFusions);
        }

        return fusions;
    }

    public final List<GeneFusion> findFusions(final List<GeneAnnotation> breakendGenes1, final List<GeneAnnotation> breakendGenes2)
    {
        final List<GeneFusion> potentialFusions = Lists.newArrayList();

        for (final GeneAnnotation startGene : breakendGenes1)
        {
            // left is upstream, right is downstream
            boolean startUpstream = isUpstream(startGene);

            for (final GeneAnnotation endGene : breakendGenes2)
            {
                boolean endUpstream = isUpstream(endGene);

                if (startUpstream == endUpstream)
                    continue;

                for (final Transcript startTrans : startGene.transcripts())
                {
                    for (final Transcript endTrans : endGene.transcripts())
                    {
                        final Transcript upstreamTrans = startUpstream ? startTrans : endTrans;
                        final Transcript downstreamTrans = !startUpstream ? startTrans : endTrans;

                        GeneFusion geneFusion = checkFusionLogic(upstreamTrans, downstreamTrans, true);

                        if(geneFusion == null)
                            continue;

                        geneFusion.setPrimarySource(
                                mKnownFusionsModel.primarySource(upstreamTrans.parent().synonyms(), downstreamTrans.parent().synonyms()));

                        potentialFusions.add(geneFusion);
                    }
                }
            }
        }

        setReportableGeneFusions(potentialFusions);

        return potentialFusions;
    }

    public static boolean validFusionTranscript(final Transcript transcript)
    {
        // check any conditions which would preclude this transcript being a part of a fusion no matter the other end
        if(transcript.postCoding())
            return false;

        boolean isUpstream = isUpstream(transcript.parent());

        if(isUpstream)
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

    public static GeneFusion checkFusionLogic(final Transcript upstreamTrans, final Transcript downstreamTrans, boolean requirePhaseMatch)
    {
        // see FV Fusions document for permitted combinations
        boolean checkExactMatch = false;

        if(!validFusionTranscript(upstreamTrans) || !validFusionTranscript(downstreamTrans))
            return null;

        if(upstreamTrans.preCoding())
        {
            if(upstreamTrans.isExonic() && !downstreamTrans.isExonic())
                return null;
            else if(downstreamTrans.isCoding())
                return null;
        }
        else if(upstreamTrans.isCoding())
        {
            if(!downstreamTrans.isCoding())
                return null;

            if(upstreamTrans.isExonic())
            {
                if(!downstreamTrans.isExonic())
                    return null;

                if(upstreamTrans.parent().id() != downstreamTrans.parent().id())
                    return null;

                // coding exon to coding exon will require phase adjustments to be exact
                checkExactMatch = true;
            }
        }
        else if(upstreamTrans.nonCoding())
        {
            if(upstreamTrans.isExonic() && !downstreamTrans.isExonic())
                return null;
            else if(downstreamTrans.isCoding())
                return null;
        }

        if (!isPotentiallyRelevantFusion(upstreamTrans, downstreamTrans))
            return null;

        boolean phaseMatched = false;

        if(!checkExactMatch)
        {
            // all fusions to downstream exons may be excluded, but for now definitely exclude those which end in the last exon
            if(downstreamTrans.isExonic() && downstreamTrans.ExonDownstream == downstreamTrans.ExonMax && !downstreamTrans.preCoding())
                return null;
        }

        if(checkExactMatch)
        {
            phaseMatched = exonToExonInPhase(upstreamTrans, true, downstreamTrans, false);

            if(phaseMatched || !requirePhaseMatch)
            {
                return new GeneFusion(upstreamTrans, downstreamTrans, phaseMatched, true);
            }
        }
        else
        {
            // just check for a phasing match
            phaseMatched = upstreamTrans.ExonUpstreamPhase == downstreamTrans.ExonDownstreamPhase;

            if(phaseMatched || !requirePhaseMatch)
            {
                return new GeneFusion(upstreamTrans, downstreamTrans, phaseMatched, true);
            }
        }

        return null;
    }

    private static boolean exonToExonInPhase(final Transcript startTrans, boolean startUpstream, final Transcript endTrans, boolean endUpstream)
    {
        // check phasing and offset since exon start or coding start
        int calcStartPhase = calcPositionPhasing(startTrans, startUpstream);

        // factor in insert sequence
        if(!startTrans.parent().insertSequence().isEmpty())
        {
            int insSeqAdjustment = (startTrans.parent().insertSequence().length() % 3);
            calcStartPhase += insSeqAdjustment;
        }

        int calcEndPhase = calcPositionPhasing(endTrans, endUpstream);

        startTrans.setExactCodingBase(calcStartPhase);
        endTrans.setExactCodingBase(calcEndPhase);

        return calcStartPhase == calcEndPhase;
    }

    private static int calcPositionPhasing(final Transcript transcript, boolean isUpstream)
    {
        // if upstream then can just use the coding bases
        // if downstream then coding bases are what's remaining
        long codingBases = transcript.calcCodingBases(isUpstream);

        int adjustedPhase = (int)(codingBases % 3);

        return adjustedPhase;
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

    private void setReportableGeneFusions(final List<GeneFusion> fusions)
    {
        GeneFusion reportableFusion = determineReportableFusion(fusions);

        if(reportableFusion == null)
            return;

        if(intragenic(reportableFusion.upstreamTrans(), reportableFusion.downstreamTrans())
        && reportableFusion.upstreamTrans().ExonUpstreamPhase == -1)
        {
            return;
        }

        reportableFusion.setReportable(true);

        // check impact on protein regions
        final Transcript downTrans = reportableFusion.downstreamTrans();
        final List<TranscriptProteinData> transProteinData = mGeneTranscriptCollection.getTranscriptProteinDataMap().get(downTrans.TransId);

        if(transProteinData == null || transProteinData.isEmpty() || !downTrans.isCoding())
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

        if (reportableFusion.getKnownFusionType() != REPORTABLE_TYPE_KNOWN)
        {
            long requiredKeptButLost = mProteinsRequiredKept.stream().filter(f -> downTrans.getProteinFeaturesLost().contains(f)).count();
            long requiredLostButKept = mProteinsRequiredLost.stream().filter(f -> downTrans.getProteinFeaturesKept().contains(f)).count();

            if (requiredKeptButLost > 0 || requiredLostButKept > 0)
                reportableFusion.setReportable(false);
        }
    }

    private void addProteinFeature(final Transcript transcript, boolean isDownstream, final String feature, int featureStart, int featureEnd)
    {
        boolean featurePreserved;

        if(isDownstream)
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

        transcript.addProteinFeature(feature, featurePreserved);
    }

    public static String TRANSCRIPT_PROTEIN_CODING = "protein_coding";
    public static String TRANSCRIPT_NONSENSE_MED_DECAY = "nonsense_mediated_decay";

    private static int MAX_UPSTREAM_DISTANCE_KNOWN = 100000;
    private static int MAX_UPSTREAM_DISTANCE_UNKNOWN = 10000;

    private GeneFusion determineReportableFusion(final List<GeneFusion> fusions)
    {
        // Select either the canonical -> canonical transcript fusion
        //  then the one with the most exons where one end is canonical
        //  then the one with the most exons combined transcript

        GeneFusion reportableFusion = null;

        // form a score by allocating 0/1 or length value to each power 10 descending
        long highestScore = 0;

        for(final GeneFusion fusion : fusions)
        {
            // first check whether a fusion is known or not - a key requirement of it being potentially reportable
            final Transcript upTrans = fusion.upstreamTrans();
            final Transcript downTrans = fusion.downstreamTrans();

            final String knownType = getKnownFusionType(upTrans, downTrans);

            if (knownType == REPORTABLE_TYPE_NONE)
                continue;

            fusion.setKnownFusionType(knownType);

            // set limits on how far upstream the breakend can be - adjusted for whether the fusions is known or not
            int maxUpstreamDistance = knownType == REPORTABLE_TYPE_KNOWN ? MAX_UPSTREAM_DISTANCE_KNOWN : MAX_UPSTREAM_DISTANCE_UNKNOWN;

            if(upTrans.getDistanceUpstream() > maxUpstreamDistance || downTrans.getDistanceUpstream() > maxUpstreamDistance)
                continue;

            if(downTrans.bioType().equals(TRANSCRIPT_NONSENSE_MED_DECAY))
                continue;

            if(downTrans.exonDistanceUp() < 0)
                continue;

            /* prioritisation rules:
            - take both canonical if possible
            - favour 3' over 5' by canoncial, protein-coding then coding bases (or exon count if not coding)
            */

            if(downTrans.isCanonical() && upTrans.isCanonical())
                return fusion;

            long transScore = 0;
            long factor = 10000000;

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


    private String getKnownFusionType(final Transcript upTrans, final Transcript downTrans)
    {
        if(mKnownFusionsModel.exactMatch(upTrans.parent().synonyms(), downTrans.parent().synonyms()))
            return REPORTABLE_TYPE_KNOWN;

        boolean intergenicPromiscuousMatch = mKnownFusionsModel.intergenicPromiscuousMatch(upTrans.parent().synonyms(), downTrans.parent().synonyms());

        boolean intragenicPromiscuousMatch = (mKnownFusionsModel.intragenicPromiscuousMatch(upTrans.parent().synonyms(), downTrans.parent().synonyms())
                && downTrans.ExonDownstream - upTrans.ExonUpstream > EXON_THRESHOLD);

        if(intergenicPromiscuousMatch || intragenicPromiscuousMatch)
        {
            boolean fivePrimeMatch = mKnownFusionsModel.fivePrimePromiscuousMatch(upTrans.parent().synonyms());
            boolean threePrimeMatch = mKnownFusionsModel.threePrimePromiscuousMatch(downTrans.parent().synonyms());

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

    public void initialiseOutputFile(final String sampleId, boolean hasMultipleSamples, String clusterInfoHeaders)
    {
        try
        {
            if(mFusionWriter == null)
            {
                String outputFilename = mOutputDir;

                if (!outputFilename.endsWith(File.separator))
                    outputFilename += File.separator;

                if(hasMultipleSamples)
                    outputFilename += "FUSIONS.csv";
                else
                    outputFilename += sampleId + "_fusions.csv";

                mFusionWriter = createBufferedWriter(outputFilename, false);

                mFusionWriter.write("SampleId,Reportable,KnownType,PrimarySource");

                if(clusterInfoHeaders.isEmpty())
                    mFusionWriter.write(",ClusterId,ClusterCount,ClusterInfo");
                else
                    mFusionWriter.write(String.format(",%s", clusterInfoHeaders));


                mFusionWriter.write(",SvIdUp,ChrUp,PosUp,OrientUp,TypeUp,PloidyUp,GeneUp,ChrBandUp,TranscriptUp,StrandUp,RegionTypeUp,CodingTypeUp");
                mFusionWriter.write(",ExonUp,PhaseUp,ExonMaxUp,DisruptiveUp,ExactBaseUp,CodingBasesUp,TotalCodingUp");
                mFusionWriter.write(",CodingStartUp,CodingEndUp,TransStartUp,TransEndUp,DistancePrevUp,CanonicalUp,BiotypeUp");

                mFusionWriter.write(",SvIdDown,ChrDown,PosDown,OrientDown,TypeDown,PloidyDown,GeneDown,ChrBandDown,TranscriptDown,StrandDown,RegionTypeDown,CodingTypeDown");
                mFusionWriter.write(",ExonDown,PhaseDown,ExonMaxDown,DisruptiveDown,ExactBaseDown,CodingBasesDown,TotalCodingDown");
                mFusionWriter.write(",CodingStartDown,CodingEndDown,TransStartDown,TransEndDown,DistancePrevDown,CanonicalDown,BiotypeDown");

                mFusionWriter.write(",ProteinsKept,ProteinsLost");
                mFusionWriter.newLine();
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing fusions: {}", e.toString());
        }
    }

    public void writeFusions(final List<GeneFusion> fusions, final String sampleId, boolean hasMultipleSamples)
    {
        initialiseOutputFile(sampleId, hasMultipleSamples, "");
        String clusterInfo = ",,";

        for (final GeneFusion fusion : fusions)
        {
            writeFusionData(fusion, sampleId, clusterInfo);
        }
    }

    public void writeFusionData(final GeneFusion fusion, final String sampleId, final String clusterInfo)
    {
        try
        {
            BufferedWriter writer = mFusionWriter;

            final Transcript upTrans = fusion.upstreamTrans();
            final Transcript downTrans = fusion.downstreamTrans();

            final GeneAnnotation startVar = upTrans.parent();
            final GeneAnnotation endVar = downTrans.parent();

            writer.write(String.format("%s,%s,%s,%s,%s",
                    sampleId, fusion.reportable(), fusion.getKnownFusionType(), fusion.primarySource(), clusterInfo));

            // write upstream SV, transcript and exon info
            writer.write(
                    String.format(",%d,%s,%d,%d,%s,%.6f",
                            startVar.id(), startVar.chromosome(), startVar.position(), startVar.orientation(),
                            startVar.type(), startVar.ploidy()));

            writer.write(
                    String.format(",%s,%s,%s,%d,%s,%s",
                            startVar.GeneName, startVar.karyotypeBand(), upTrans.StableId,
                            startVar.Strand, upTrans.regionType(), upTrans.codingType()));

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
                            endVar.id(), endVar.chromosome(), endVar.position(), endVar.orientation(),
                            endVar.type(), endVar.ploidy()));

            writer.write(
                    String.format(",%s,%s,%s,%d,%s,%s",
                            endVar.GeneName, endVar.karyotypeBand(), downTrans.StableId,
                            endVar.Strand, downTrans.regionType(), downTrans.codingType()));

            writer.write(
                    String.format(",%d,%d,%d,%s",
                            downTrans.ExonDownstream, downTrans.ExonDownstreamPhase, downTrans.ExonMax, downTrans.isDisruptive()));

            writer.write(
                    String.format(",%d,%d,%d,%d,%d,%d,%d,%d,%s,%s,%s,%s",
                            downTrans.exactCodingBase(), downTrans.calcCodingBases(false), downTrans.totalCodingBases(),
                            downTrans.codingStart(), downTrans.codingEnd(), downTrans.TranscriptStart, downTrans.TranscriptEnd,
                            downTrans.exonDistanceUp(), downTrans.isCanonical(), downTrans.bioType(),
                            downTrans.getProteinFeaturesKept(), downTrans.getProteinFeaturesLost()));

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
        final List<TranscriptExonData> exonDataList = mGeneTranscriptCollection.getTranscriptExons(trans.parent().StableId, trans.StableId);

        if (exonDataList == null || exonDataList.isEmpty())
            return false;

        int strand = trans.parent().Strand;

        // first find the matching exon boundary for this RNA fusion boundary
        for (int i = 0; i < exonDataList.size(); ++i)
        {
            final TranscriptExonData exonData = exonDataList.get(i);
            final TranscriptExonData prevExonData = i > 0 ? exonDataList.get(i - 1) : null;
            final TranscriptExonData nextExonData = i < exonDataList.size() - 1 ? exonDataList.get(i + 1) : null;

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
