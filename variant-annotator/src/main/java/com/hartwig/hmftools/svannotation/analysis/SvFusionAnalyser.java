package com.hartwig.hmftools.svannotation.analysis;

import static com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalyzer.getIntronicTranscripts;
import static com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalyzer.intronicDisruptionOnSameTranscript;
import static com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalyzer.isUpstream;

import java.util.Comparator;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;

import org.apache.commons.cli.Options;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class SvFusionAnalyser
{

    public static final String FUSION_PAIRS_CSV = "fusion_pairs_csv";
    public static final String PROMISCUOUS_FIVE_CSV = "promiscuous_five_csv";
    public static final String PROMISCUOUS_THREE_CSV = "promiscuous_three_csv";

    private static final int EXON_THRESHOLD = 1;

    private final KnownFusionsModel mKnownFusionsModel;

    private static final Logger LOGGER = LogManager.getLogger(SvFusionAnalyser.class);

    public SvFusionAnalyser(final KnownFusionsModel knownFusionsModel)
    {
        mKnownFusionsModel = knownFusionsModel;
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
        final List<List<Pair<Transcript, Transcript>>> fusionTranscriptPairs = Lists.newArrayList();

        final List<Pair<Transcript, Transcript>> potentialFusions = Lists.newArrayList();

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
                    if(startTrans.isPromoter())
                        continue;

                    for (final Transcript endTrans : endGene.transcripts())
                    {
                        if(endTrans.isPromoter())
                            continue;

                        if (!isPotentiallyRelevantFusion(startTrans, endTrans))
                            continue;

                        if(startTrans.isIntronic())
                        {
                            // Intron -> Intron and Intron -> Exon are both fine if phased
                            if (startUpstream && startTrans.exonUpstreamPhase() == endTrans.exonDownstreamPhase())
                            {
                                addFusion(potentialFusions, startTrans, endTrans);
                            }
                            else if (!startUpstream && endTrans.exonUpstreamPhase() == startTrans.exonDownstreamPhase())
                            {
                                addFusion(potentialFusions, endTrans, startTrans);
                            }
                        }
                        else if(startTrans.isExonic())
                        {
                            if(endTrans.isExonic())
                            {
                                if(exonToExonInPhase(startTrans, startUpstream, endTrans, endUpstream))
                                {
                                    addFusion(potentialFusions, startTrans, endTrans);
                                }
                                else if(exonToExonInPhase(endTrans, endUpstream, startTrans, startUpstream))
                                {
                                    addFusion(potentialFusions, endTrans, startTrans);
                                }
                            }

                            // Exon -> Intron is invalid

                        }
                        else
                        {
                            // UTR region fusions not handled yet
                        }
                    }
                }
            }
        }

        fusionTranscriptPairs.add(potentialFusions);

        List<GeneFusion> fusions = toReportableGeneFusions(fusionTranscriptPairs);

        return fusions;
    }

    private static boolean exonToExonInPhase(final Transcript startTrans, boolean startUpstream, final Transcript endTrans, boolean endUpstream)
    {
        // check phasing and offset since exon start or coding start
        long calcStartPhase = calcPositionPhasing(startTrans, startUpstream);
        long calcEndPhase = calcPositionPhasing(endTrans, endUpstream);

        return calcStartPhase == calcEndPhase;
    }

    private static long calcPositionPhasing(final Transcript transcript, boolean isUpstream)
    {
        // if upstream then can just use the coding bases
        // if downstream then coding bases are what's remaing
        long codingBases = isUpstream ? transcript.codingBases() : transcript.totalCodingBases() - transcript.codingBases();

        int adjustedPhase = (int)(codingBases % 3);

        return adjustedPhase;
    }

    private static int calcPositionPhasing_v2(final Transcript transcript, boolean isUpstream)
    {
        // if the exon is completely translated, then determine bases since start of exon
        long position = transcript.parent().position();

        long exonOffset = 0;

        if(transcript.codingStart() != null && position <= transcript.codingStart())
        {
            exonOffset = -1;
        }
        else if(transcript.codingStart() != null && position > transcript.codingStart())
        {
            exonOffset = position - transcript.codingStart();
        }
        if(transcript.codingEnd() != null && position >= transcript.codingEnd())
        {
            exonOffset = -1;
        }
        /*
        else if(transcript.codingEnd() != null && position < transcript.codingEnd())
        {
            exonOffset = position - transcript.exonStart();
        }
        */

        int phasing = isUpstream ? transcript.exonUpstreamPhase() : transcript.exonDownstreamPhase();

        int combinedPhasing = (int)(phasing + exonOffset);
        int adjustedPhase = (combinedPhasing % 3);

        return adjustedPhase;
    }

    private void addFusion(List<Pair<Transcript, Transcript>> fusions, final Transcript startTrans, final Transcript endTrans)
    {
        //LOGGER.debug("adding fusion between start SV({}) trans({}) and end SV({}) trans({})",
        //        startTrans.parent().id(), startTrans.toString(), endTrans.parent().id(), endTrans.toString());

        fusions.add(Pair.of(startTrans, endTrans));
    }

    private static boolean isPotentiallyRelevantFusion(final Transcript t1, final Transcript t2)
    {
        final boolean sameGene = t1.geneName().equals(t2.geneName());

        if (sameGene)
        {
            // skip fusions between different transcripts in the same gene,
            if (!t1.transcriptId().equals(t2.transcriptId()))
                return false;

            // skip fusions within the same intron
            if (intronicDisruptionOnSameTranscript(t1, t2))
                return false;
        }

        return true;
    }

    @NotNull
    private List<GeneFusion> toReportableGeneFusions(@NotNull List<List<Pair<Transcript, Transcript>>> fusionTranscriptPairs)
    {
        final List<GeneFusion> result = Lists.newArrayList();

        for (final List<Pair<Transcript, Transcript>> fusions : fusionTranscriptPairs)
        {
            Optional<Pair<Transcript, Transcript>> reportableFusion = determineReportableFusion(fusions);

            for (final Pair<Transcript, Transcript> fusion : fusions)
            {
                final Transcript upstream = fusion.getLeft();
                final Transcript downstream = fusion.getRight();
                final boolean matchesKnownFusion = transcriptsMatchKnownFusion(mKnownFusionsModel, upstream, downstream);

                Boolean isPostCodingUpstream = postCoding(upstream);
                Boolean isPostCodingDownstream = postCoding(downstream);

                final boolean reportable = reportableFusion.isPresent() && reportableFusion.get() == fusion && matchesKnownFusion && (
                        (isPostCodingDownstream != null && !isPostCodingDownstream) && (isPostCodingUpstream == null
                                || !isPostCodingUpstream) && !(intragenic(upstream, downstream) && upstream.exonUpstreamPhase() == -1));

                final GeneFusion geneFusion = ImmutableGeneFusion.builder()
                        .reportable(reportable)
                        .upstreamLinkedAnnotation(upstream)
                        .downstreamLinkedAnnotation(downstream)
                        .primarySource(mKnownFusionsModel.primarySource(upstream.parent().synonyms(), downstream.parent().synonyms()))
                        .build();

                result.add(geneFusion);
            }
        }
        return result;
    }

    @NotNull
    private static Optional<Pair<Transcript, Transcript>> determineReportableFusion(@NotNull List<Pair<Transcript, Transcript>> fusions)
    {
        // Select either the canonical -> canonical transcript fusion
        //  then the one with the most exons where one end is canonical
        //  then the one with the most exons combined transcript

        Optional<Pair<Transcript, Transcript>> reportableFusion =
                fusions.stream().filter(pair -> pair.getLeft().isCanonical() && pair.getRight().isCanonical()).findFirst();

        if (!reportableFusion.isPresent())
        {
            reportableFusion = fusions.stream()
                    .filter(pair -> pair.getLeft().isCanonical() || pair.getRight().isCanonical())
                    .sorted(Comparator.comparingInt(a -> a.getLeft().exonMax() + a.getRight().exonMax()))
                    .reduce((a, b) -> b);
        }

        if (!reportableFusion.isPresent())
        {
            reportableFusion = fusions.stream()
                    .sorted(Comparator.comparingInt(a -> a.getLeft().exonMax() + a.getRight().exonMax()))
                    .reduce((a, b) -> b);
        }

        return reportableFusion;
    }


    private static boolean transcriptsMatchKnownFusion(final KnownFusionsModel fusionsModel, final Transcript five, final Transcript three)
    {
        if(fusionsModel.exactMatch(five.parent().synonyms(), three.parent().synonyms()))
            return true;

        if(fusionsModel.intergenicPromiscuousMatch(five.parent().synonyms(), three.parent().synonyms()))
            return true;

        if(fusionsModel.intragenicPromiscuousMatch(five.parent().synonyms(), three.parent().synonyms())
        && three.exonDownstream() - five.exonUpstream() > EXON_THRESHOLD)
        {
            return true;
        }

        return false;
    }

    private static Boolean postCoding(@NotNull final Transcript transcript)
    {
        Long codingStart = transcript.codingStart();
        Long codingEnd = transcript.codingEnd();

        if (codingStart == null || codingEnd == null)
            return null;

        final int strand = transcript.parent().strand();
        final long position = transcript.parent().position();

        return (strand == 1 && position > codingEnd) || (strand == -1 && position < codingStart);
    }

    private static boolean intragenic(final Transcript upstream, final Transcript downstream)
    {
        return upstream.parent().synonyms().stream().anyMatch(downstream.parent().synonyms()::contains);
    }

}
