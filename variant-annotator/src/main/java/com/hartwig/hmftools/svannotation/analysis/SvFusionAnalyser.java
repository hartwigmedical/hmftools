package com.hartwig.hmftools.svannotation.analysis;

import static com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalyzer.isUpstream;

import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;

import org.apache.commons.cli.Options;
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

    boolean mIncludePossibles;

    private static final Logger LOGGER = LogManager.getLogger(SvFusionAnalyser.class);

    public SvFusionAnalyser(final KnownFusionsModel knownFusionsModel)
    {
        mKnownFusionsModel = knownFusionsModel;
        mIncludePossibles = false;
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(FUSION_PAIRS_CSV, true, "Path towards a CSV containing white-listed gene fusion pairs.");
        options.addOption(PROMISCUOUS_FIVE_CSV, true, "Path towards a CSV containing white-listed promiscuous 5' genes.");
        options.addOption(PROMISCUOUS_THREE_CSV, true, "Path towards a CSV containing white-listed promiscuous 3' genes.");
    }

    public void setIncludePossibles(boolean toggle) { mIncludePossibles = toggle; }

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

                // see FV Fusions document for permitted combinations
                for (final Transcript startTrans : startGene.transcripts())
                {
                    for (final Transcript endTrans : endGene.transcripts())
                    {
                        final Transcript upstreamTrans = startUpstream ? startTrans : endTrans;
                        final Transcript downstreamTrans = !startUpstream ? startTrans : endTrans;

                        boolean checkExactMatch = false;

                        if(upstreamTrans.postCoding() || downstreamTrans.postCoding() || downstreamTrans.nonCoding())
                            continue;

                        if(upstreamTrans.isPromoter())
                            continue;

                        if(downstreamTrans.isPrePromotor())
                            continue;

                        if(upstreamTrans.preCoding())
                        {
                            if(upstreamTrans.isExonic() && !downstreamTrans.isExonic())
                                continue;
                            else if(downstreamTrans.isCoding())
                                continue;

                            // phasing match
                        }
                        else if(upstreamTrans.isCoding())
                        {
                            if(!downstreamTrans.isCoding())
                                continue;

                            if(upstreamTrans.isExonic())
                            {
                                if(!downstreamTrans.isExonic())
                                    continue;

                                checkExactMatch = true;
                            }

                            // phasing match
                        }
                        else if(upstreamTrans.nonCoding())
                        {
                            if(upstreamTrans.isExonic() && !downstreamTrans.isExonic())
                                continue;
                            else if(downstreamTrans.isCoding())
                                continue;

                            // phasing match
                        }

                        /*
                        // DEBUG
                        if(upstreamTrans.transcriptId().equals("ENST00000312970") && downstreamTrans.transcriptId().equals("ENST00000438429"))
                        {
                            LOGGER.debug("trans match");
                        }
                        */

                        if (!isPotentiallyRelevantFusion(upstreamTrans, downstreamTrans))
                            continue;

                        if(!checkExactMatch)
                        {
                            // all fusions to downstream exons may be excluded, but for now definitely exclude those which end in the last exon
                            if(downstreamTrans.isExonic() && downstreamTrans.exonDownstream() == downstreamTrans.exonMax() && ! downstreamTrans.preCoding())
                                continue;
                        }

                        if(checkExactMatch)
                        {
                            if(exonToExonInPhase(upstreamTrans, true, downstreamTrans, false))
                            {
                                addFusion(potentialFusions, upstreamTrans, downstreamTrans, true);
                            }
                        }
                        else
                        {
                            // just check for a phasing match
                            if (upstreamTrans.exonUpstreamPhase() == downstreamTrans.exonDownstreamPhase())
                            {
                                addFusion(potentialFusions, upstreamTrans, downstreamTrans, true);
                            }
                        }

                        if(mIncludePossibles && transcriptsMatchKnownFusion(upstreamTrans, downstreamTrans))
                        {
                            addFusion(potentialFusions, upstreamTrans, downstreamTrans, false);
                        }
                    }
                }
            }
        }

        setReportableGeneFusions(potentialFusions);

        return potentialFusions;
    }

    private static boolean exonToExonInPhase(final Transcript startTrans, boolean startUpstream, final Transcript endTrans, boolean endUpstream)
    {
        // check phasing and offset since exon start or coding start
        int calcStartPhase = calcPositionPhasing(startTrans, startUpstream);
        int calcEndPhase = calcPositionPhasing(endTrans, endUpstream);

        startTrans.setExactCodingBase(calcStartPhase);
        endTrans.setExactCodingBase(calcEndPhase);

        return calcStartPhase == calcEndPhase;
    }

    private static int calcPositionPhasing(final Transcript transcript, boolean isUpstream)
    {
        // if upstream then can just use the coding bases
        // if downstream then coding bases are what's remaing
        long codingBases = isUpstream ? transcript.codingBases() : transcript.totalCodingBases() - transcript.codingBases();

        int adjustedPhase = (int)(codingBases % 3);

        return adjustedPhase;
    }

    private void addFusion(List<GeneFusion> fusions, final Transcript upstreamTrans, final Transcript downstreamTrans, boolean phaseMatched)
    {
        //LOGGER.debug("adding fusion between start SV({}) trans({}) and end SV({}) trans({})",
        //        startTrans.parent().id(), startTrans.toString(), endTrans.parent().id(), endTrans.toString());

        fusions.add(new GeneFusion(
                upstreamTrans,
                downstreamTrans,
                mKnownFusionsModel.primarySource(upstreamTrans.parent().synonyms(), downstreamTrans.parent().synonyms()),
                false,
                phaseMatched));
    }

    private static boolean isPotentiallyRelevantFusion(final Transcript t1, final Transcript t2)
    {
        if(!t1.geneName().equals(t2.geneName()))
            return true;

        // skip fusions between different transcripts in the same gene,
        if (!t1.transcriptId().equals(t2.transcriptId()))
            return false;

        if(t1.nonCoding())
            return false;

        // skip fusions within the same intron
        if(t1.isIntronic() && t2.isIntronic() && t1.exonUpstream() == t2.exonUpstream())
            return false;

        return true;
    }

    private void setReportableGeneFusions(final List<GeneFusion> fusions)
    {
        final List<GeneFusion> result = Lists.newArrayList();

        Optional<GeneFusion> reportableFusion = determineReportableFusion(fusions);

        if(!reportableFusion.isPresent() )
            return;

        for (final GeneFusion fusion : fusions)
        {
            if(!fusion.isPhaseMatch())
                continue;

            if(reportableFusion.get() == fusion)
            {
                boolean intragenicOk = !(intragenic(fusion.upstreamTrans(), fusion.downstreamTrans()) && fusion.upstreamTrans().exonUpstreamPhase() == -1);

                if(intragenicOk)
                    fusion.setReportable(true);

                // old rule was checking (downstream: not non-coding and not downstream from coding)
                // AND (upstream: non-coding or not downstream of coding)
                // AND NOT (intragenic and upstream phase -1)
                /*
                boolean oldReportable = isPostCodingDownstream != null && !isPostCodingDownstream
                        && (isPostCodingUpstream == null || !isPostCodingUpstream)
                        && !(intragenic(upstream, downstream) && upstream.exonUpstreamPhase() == -1);
                */
            }
        }

    }

    @NotNull
    private Optional<GeneFusion> determineReportableFusion(final List<GeneFusion> fusions)
    {
        // Select either the canonical -> canonical transcript fusion
        //  then the one with the most exons where one end is canonical
        //  then the one with the most exons combined transcript

        List<GeneFusion> knownFusions = fusions.stream()
                .filter(f -> transcriptsMatchKnownFusion(f.upstreamTrans(), f.downstreamTrans()))
                .collect(Collectors.toList());

        Optional<GeneFusion> reportableFusion =
                knownFusions.stream()
                        .filter(f -> f.upstreamTrans().isCanonical() && f.downstreamTrans().isCanonical()).findFirst();

        if (!reportableFusion.isPresent())
        {
            reportableFusion = knownFusions.stream()
                    .filter(f -> f.upstreamTrans().isCanonical() || f.downstreamTrans().isCanonical())
                    .sorted(Comparator.comparingInt(a -> a.upstreamTrans().exonMax() + a.downstreamTrans().exonMax()))
                    .reduce((a, b) -> b);
        }

        if (!reportableFusion.isPresent())
        {
            reportableFusion = knownFusions.stream()
                    .sorted(Comparator.comparingInt(a -> a.upstreamTrans().exonMax() + a.downstreamTrans().exonMax()))
                    .reduce((a, b) -> b);
        }

        return reportableFusion;
    }

    private boolean transcriptsMatchKnownFusion(final Transcript upTrans, final Transcript downTrans)
    {
        if(mKnownFusionsModel.exactMatch(upTrans.parent().synonyms(), downTrans.parent().synonyms()))
            return true;

        if(mKnownFusionsModel.intergenicPromiscuousMatch(upTrans.parent().synonyms(), downTrans.parent().synonyms()))
            return true;

        if(mKnownFusionsModel.intragenicPromiscuousMatch(upTrans.parent().synonyms(), downTrans.parent().synonyms())
        && downTrans.exonDownstream() - upTrans.exonUpstream() > EXON_THRESHOLD)
        {
            return true;
        }

        return false;
    }

    private static boolean intragenic(final Transcript upstream, final Transcript downstream)
    {
        return upstream.parent().synonyms().stream().anyMatch(downstream.parent().synonyms()::contains);
    }

}
