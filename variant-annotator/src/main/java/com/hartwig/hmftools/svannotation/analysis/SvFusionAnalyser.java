package com.hartwig.hmftools.svannotation.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData.GENE_PHASING_REGION_5P_UTR;
import static com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData.GENE_PHASING_REGION_CODING_0;
import static com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData.GENE_PHASING_REGION_CODING_1;
import static com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData.GENE_PHASING_REGION_CODING_2;
import static com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData.GENE_PHASING_REGION_MAX;
import static com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData.phaseToRegion;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation.isDownstream;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation.isUpstream;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.EXON_RANK_MAX;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.EXON_RANK_MIN;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.nextTranscriptExons;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptExonData;
import com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection;

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

    private Map<String, List<RnaFusionData>> mSampleRnaData;
    private SvGeneTranscriptCollection mGeneTranscriptCollection;

    boolean mIncludePossibles;

    private BufferedWriter mFusionWriter;
    private BufferedWriter mRnaWriter;

    private static final Logger LOGGER = LogManager.getLogger(SvFusionAnalyser.class);

    public SvFusionAnalyser(final KnownFusionsModel knownFusionsModel, final SvGeneTranscriptCollection geneTranscriptCollection)
    {
        mKnownFusionsModel = knownFusionsModel;
        mIncludePossibles = false;
        mFusionWriter = null;
        mRnaWriter = null;
        mSampleRnaData = new HashMap();
        mGeneTranscriptCollection = geneTranscriptCollection;
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

    private static int SPECIFIC_VAR_ID = -1;
    // private static int SPECIFIC_VAR_ID = 5232794;

    public final List<GeneFusion> findFusions(final List<GeneAnnotation> breakendGenes1, final List<GeneAnnotation> breakendGenes2)
    {
        final List<GeneFusion> potentialFusions = Lists.newArrayList();

        for (final GeneAnnotation startGene : breakendGenes1)
        {
            if(startGene.id() == SPECIFIC_VAR_ID)
            {
                LOGGER.debug("specific var({})", startGene.id());
            }

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

                        /*
                        // DEBUG
                        if(upstreamTrans.transcriptId().equals("ENST00000328159") && downstreamTrans.transcriptId().equals("ENST00000392403"))
                        {
                            LOGGER.debug("trans match");
                        }
                         */

                        boolean checkExactMatch = false;

                        if(upstreamTrans.postCoding() || downstreamTrans.postCoding() || downstreamTrans.nonCoding())
                            continue;

                        if(upstreamTrans.isPromoter())
                            continue;

                        if(downstreamTrans.isPrePromotor())
                            continue;

                        if(downstreamTrans.exonMax() == 1)
                            continue;

                        if(!upstreamTrans.isDisruptive())
                            continue;;

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

                                // skip chained fusions from exon-to-exon?
                                if(upstreamTrans.parent().id() != downstreamTrans.parent().id())
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

        // factor in insert sequence
        if(!startTrans.parent().insertSequence().isEmpty())
        {
            int insSeqAdjustment = (int)(startTrans.parent().insertSequence().length() % 3);
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
        GeneFusion reportableFusion = determineReportableFusion(fusions);

        if(reportableFusion == null)
            return;

        for (final GeneFusion fusion : fusions)
        {
            if(!fusion.isPhaseMatch())
                continue;

            if(reportableFusion == fusion)
            {
                boolean intragenicOk = !(intragenic(fusion.upstreamTrans(), fusion.downstreamTrans()) && fusion.upstreamTrans().exonUpstreamPhase() == -1);

                if(intragenicOk)
                    fusion.setReportable(true);
            }
        }
    }

    public static String TRANSCRIPT_PROTEIN_CODING = "protein_coding";
    public static String TRANSCRIPT_NONSENSE_MED_DECAY = "nonsense_mediated_decay";

    private GeneFusion determineReportableFusion(final List<GeneFusion> fusions)
    {
        // Select either the canonical -> canonical transcript fusion
        //  then the one with the most exons where one end is canonical
        //  then the one with the most exons combined transcript

        List<GeneFusion> knownFusions = fusions.stream()
                .filter(GeneFusion::isPhaseMatch)
                .filter(f -> transcriptsMatchKnownFusion(f.upstreamTrans(), f.downstreamTrans()))
                .filter(f -> !f.downstreamTrans().bioType().equals(TRANSCRIPT_NONSENSE_MED_DECAY))
                .filter(f -> f.downstreamTrans().exonDistanceUp() >= 0)
                .collect(Collectors.toList());

        if(knownFusions.isEmpty())
            return null;

        /* new prioritisation rules:
        - take both canonical if possible
        - favour 3' over 5' by canoncial, protein-coding then coding bases (or exon count if not coding)
        */

        GeneFusion reportableFusion = null;

        // form a score by allocating 0/1 or length value to each power 10 descending
        long highestScore = 0;

        for(final GeneFusion fusion : knownFusions)
        {
            if(fusion.downstreamTrans().isCanonical() && fusion.upstreamTrans().isCanonical())
                return fusion;

            long transScore = 0;
            long factor = 10000000;

            if(fusion.downstreamTrans().isCanonical())
                transScore += factor;

            factor /= 10;

            if(fusion.downstreamTrans().bioType().equals(TRANSCRIPT_PROTEIN_CODING))
                transScore += factor;

            factor /= 100;

            long length = fusion.downstreamTrans().isCoding() ? fusion.downstreamTrans().codingBases() : fusion.downstreamTrans().exonMax();

            // will be a range between 1-99 * current factor
            length = min(round(length/10), 99);
            transScore += length * factor;

            factor /= 10;

            if(fusion.upstreamTrans().isCanonical())
                transScore += factor;

            factor /= 10;

            if(fusion.upstreamTrans().bioType().equals(TRANSCRIPT_PROTEIN_CODING))
                transScore += factor;

            factor /= 100;

            length = fusion.upstreamTrans().isCoding() ? fusion.upstreamTrans().codingBases() : fusion.upstreamTrans().exonMax();
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

    public void writeFusions(final List<GeneFusion> fusions, final String outputDir, final String sampleId,  final String clusterInfo)
    {
        if(fusions.isEmpty())
            return;

        try
        {
            if(mFusionWriter == null)
            {
                String outputFilename = outputDir;

                if (!outputFilename.endsWith("/"))
                    outputFilename += File.separator;

                outputFilename += "FUSIONS.csv";

                mFusionWriter = createBufferedWriter(outputFilename, false);

                mFusionWriter.write("SampleId,Reportable,PrimarySource,ClusterId,ClusterCount,ResolvedType,PhaseMatched");

                mFusionWriter.write(",SvIdUp,ChrUp,PosUp,OrientUp,TypeStart,GeneUp,TranscriptUp,StrandUp,RegionTypeUp,CodingTypeUp");
                mFusionWriter.write(",ExonUp,PhaseUp,ExonMaxUp,DisruptiveUp,ExactBaseUp,CodingBasesUp,TotalCodingUp");
                mFusionWriter.write(",CodingStartUp,CodingEndUp,TransStartUp,DistancePrevUp,BiotypeUp");

                mFusionWriter.write(",SvIdDown,ChrDown,PosDown,OrientDown,TypeDown,GeneDown,TranscriptDown,StrandDown,RegionTypeDown,CodingTypeDown");
                mFusionWriter.write(",ExonDown,PhaseDown,ExonMaxDown,DisruptiveDown,ExactBaseDown,CodingBasesDown,TotalCodingDown");
                mFusionWriter.write(",CodingStartDown,CodingEndDown,TransStartDown,DistancePrevDown,BiotypeDown");
                mFusionWriter.newLine();
            }

            BufferedWriter writer = mFusionWriter;

            for(final GeneFusion fusion : fusions)
            {
                final Transcript startTrans = fusion.upstreamTrans();
                final Transcript endTrans = fusion.downstreamTrans();

                final GeneAnnotation startVar = startTrans.parent();
                final GeneAnnotation endVar = endTrans.parent();

                writer.write(String.format("%s,%s,%s,%s,%s",
                        sampleId, fusion.reportable(), fusion.primarySource(), clusterInfo, fusion.isPhaseMatch()));

                // write upstream SV, transcript and exon info
                writer.write(
                        String.format(",%d,%s,%d,%d,%s",
                                startVar.id(), startVar.chromosome(), startVar.position(), startVar.orientation(), startVar.type()));

                writer.write(
                        String.format(",%s,%s,%d,%s,%s",
                                startTrans.parent().geneName(), startTrans.transcriptId(),
                                startTrans.parent().strand(), startTrans.regionType(), startTrans.codingType()));

                writer.write(
                        String.format(",%d,%d,%d,%s",
                                startTrans.exonUpstream(), startTrans.exonUpstreamPhase(), startTrans.exonMax(), startTrans.isDisruptive()));
                writer.write(
                        String.format(",%d,%d,%d,%d,%d,%d,%d,%s",
                                startTrans.exactCodingBase(), startTrans.codingBases(), startTrans.totalCodingBases(),
                                startTrans.codingStart(), startTrans.codingEnd(),
                                startTrans.transcriptStart(), startTrans.exonDistanceUp(), startTrans.bioType()));

                writer.write(
                        String.format(",%d,%s,%d,%d,%s",
                                endVar.id(), endVar.chromosome(), endVar.position(), endVar.orientation(), endVar.type()));

                writer.write(
                        String.format(",%s,%s,%d,%s,%s",
                                endTrans.parent().geneName(), endTrans.transcriptId(),
                                endTrans.parent().strand(), endTrans.regionType(), endTrans.codingType()));

                writer.write(
                        String.format(",%d,%d,%d,%s",
                                endTrans.exonDownstream(), endTrans.exonDownstreamPhase(), endTrans.exonMax(), endTrans.isDisruptive()));

                writer.write(
                        String.format(",%d,%d,%d,%d,%d,%d,%d,%s",
                                endTrans.exactCodingBase(), endTrans.codingBases(), endTrans.totalCodingBases(),
                                endTrans.codingStart(), endTrans.codingEnd(), endTrans.transcriptStart(), endTrans.exonDistanceUp(), endTrans.bioType()));

                writer.newLine();
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing fusions: {}", e.toString());
        }
    }

    public final Map<String, List<RnaFusionData>> getSampleRnaData() { return mSampleRnaData; }
    public final List<RnaFusionData> getSampleRnaData(final String sampleId) { return mSampleRnaData.get(sampleId); }

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

                RnaFusionData rnaData = new RnaFusionData(
                        items[COL_NAME],
                        items[COL_GENE_UP], items[COL_GENE_DOWN], items[COL_CHR_UP], items[COL_CHR_DOWN],
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

    public void writeRnaMatchData(final String sampleId, final String outputDir, final List<GeneFusion> fusions,
            final List<StructuralVariantAnnotation> svAnnotations)
    {
        final List<RnaFusionData> rnaFusionList = mSampleRnaData.get(sampleId);

        if(rnaFusionList == null || rnaFusionList.isEmpty())
            return;

        try
        {
            if(mRnaWriter == null)
            {
                String outputFilename = outputDir;

                if (!outputFilename.endsWith("/"))
                    outputFilename += File.separator;

                outputFilename += "RNA_MATCH_DATA.csv";

                mRnaWriter = createBufferedWriter(outputFilename, false);

                mRnaWriter.write("SampleId,FusionName,GeneUp,GeneDown,SvFusionMatch,Reportable");
                mRnaWriter.write(",ExonMinRankUp,ExonMaxRankUp,ExonMinRankDown,ExonMaxRankDown");

                mRnaWriter.write(",SvIdUp,ChrUp,PosUp,RnaPosUp,OrientUp,StrandUp,TypeStart,TranscriptUp,RegionTypeUp,CodingTypeUp");
                mRnaWriter.write(",ExonUp,PhaseUp,DisruptiveUp,TransStartUp,DistancePrevUp");

                mRnaWriter.write(",SvIdDown,ChrDown,PosDown,RnaPosDown,OrientDown,StrandDown,TypeStart,TranscriptDown,RegionTypeDown,CodingTypeDown");
                mRnaWriter.write(",ExonDown,PhaseDown,DisruptiveDown,TransStartDown,DistancePrevDown");

                mRnaWriter.write(",JunctionReadCount,SpanningFragCount,SpliceType");

                mRnaWriter.newLine();
            }


            BufferedWriter writer = mRnaWriter;

            for(final RnaFusionData rnaFusion : rnaFusionList)
            {
                matchRnaFusion(rnaFusion, fusions, svAnnotations);

                final GeneFusion fusionMatch = rnaFusion.getMatchedFusion();

                writer.write(
                        String.format("%s,%s,%s,%s,%s,%s",
                                sampleId, rnaFusion.Name, rnaFusion.GeneUp, rnaFusion.GeneDown,
                                fusionMatch != null, fusionMatch != null ? fusionMatch.reportable() : "false"));

                writer.write(String.format(",%d,%d,%d,%d",
                        rnaFusion.exonMinRankUp(), rnaFusion.exonMaxRankUp(), rnaFusion.exonMinRankDown(), rnaFusion.exonMaxRankDown()));

                final Transcript transUp = rnaFusion.getTransUp();

                if(transUp != null)
                {
                    writer.write(
                            String.format(",%d,%s,%d,%d,%d,%d,%s",
                                    transUp.parent().id(), transUp.parent().chromosome(), transUp.parent().position(), rnaFusion.PositionUp,
                                    transUp.parent().orientation(), transUp.parent().strand(), transUp.parent().type()));

                    writer.write(
                            String.format(",%s,%s,%s,%d,%d,%s,%d,%d",
                                    transUp.transcriptId(), transUp.regionType(), transUp.codingType(),
                                    transUp.exonUpstream(), transUp.exonUpstreamPhase(), transUp.isDisruptive(),
                                    transUp.transcriptStart(), transUp.exonDistanceUp()));
                }
                else
                {
                    writer.write(String.format(",,%s,,%d,,%d,", rnaFusion.ChrUp, rnaFusion.PositionUp, rnaFusion.StrandUp));
                    writer.write(",,,,,,,,");
                }

                final Transcript transDown = rnaFusion.getTransDown();

                if(transDown != null)
                {
                    writer.write(
                            String.format(",%d,%s,%d,%d,%d,%d,%s",
                                    transDown.parent().id(), transDown.parent().chromosome(), transDown.parent().position(), rnaFusion.PositionDown,
                                    transDown.parent().orientation(), transDown.parent().strand(), transDown.parent().type()));

                    writer.write(
                            String.format(",%s,%s,%s,%d,%d,%s,%d,%d",
                                    transDown.transcriptId(), transDown.regionType(), transDown.codingType(),
                                    transDown.exonDownstream(), transDown.exonDownstreamPhase(), transDown.isDisruptive(),
                                    transDown.transcriptStart(), transDown.exonDistanceUp()));
                }
                else
                {
                    writer.write(String.format(",,%s,,%d,,%d,", rnaFusion.ChrDown, rnaFusion.PositionDown, rnaFusion.StrandDown));
                    writer.write(",,,,,,,,");
                }

                writer.write(String.format(",%d,%d,%s", rnaFusion.JunctionReadCount, rnaFusion.SpanningFragCount, rnaFusion.SpliceType));

                writer.newLine();

            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing RNA match data: {}", e.toString());
        }

    }

    private void matchRnaFusion(final RnaFusionData rnaFusion, final List<GeneFusion> fusions, final List<StructuralVariantAnnotation> annotations)
    {
        // find the min and max exon rank for every transcript matching the RNA positions
        setRnaFusionData(rnaFusion);

        // first check for a fusion found
        for(final GeneFusion fusion : fusions)
        {
            final Transcript upTrans = fusion.upstreamTrans();
            final Transcript downTrans = fusion.downstreamTrans();

            if(!upTrans.geneName().equals(rnaFusion.GeneUp)  || !downTrans.geneName().equals(rnaFusion.GeneDown))
                continue;

            // breakends must be downstream of the upstream RNA exons, and upstream of the downstream RNA exon
            if(!withinPositionRange(upTrans, rnaFusion.PositionUp, upTrans.parent().strand() == 1))
                continue;

            if(!withinPositionRange(downTrans, rnaFusion.PositionDown, downTrans.parent().strand() != 1))
                continue;

            rnaFusion.setMatchedFusion(fusion);
            break;
        }

        if(rnaFusion.getMatchedFusion() == null)
        {
            Transcript transUp = null;
            Transcript transDown = null;
            boolean singleSvFound = false;

            // check breakend annotations for any match on the gene - favour a single SV which explains up and down genes
            for(final StructuralVariantAnnotation svAnnotation : annotations)
            {
                for (final GeneAnnotation breakendGene : svAnnotation.annotations())
                {
                    Transcript svTransUp = null;
                    Transcript svTransDown = null;

                    for(final Transcript transcript : breakendGene.transcripts())
                    {
                        if(isUpstream(breakendGene) && transcript.geneName().equals(rnaFusion.GeneUp)
                        && withinPositionRange(transcript, rnaFusion.PositionUp, transcript.parent().strand() == 1))
                        {
                            transUp = transcript;
                            svTransUp = transcript;
                        }

                        if(isDownstream(breakendGene) && transcript.geneName().equals(rnaFusion.GeneDown)
                        && withinPositionRange(transcript, rnaFusion.PositionDown, transcript.parent().strand() != 1))
                        {
                            transDown= transcript;
                            svTransDown = transcript;
                        }
                    }

                    if(svTransDown != null && svTransUp != null)
                    {
                        rnaFusion.setBreakends(svTransUp, svTransDown);
                        singleSvFound = true;
                    }
                }

                if(singleSvFound)
                    break;
            }

            if(!singleSvFound)
            {
                rnaFusion.setBreakends(transUp, transDown);
            }
        }
    }

    private boolean withinPositionRange(final Transcript trans, final long rnaPosition, boolean requireHigherBreakendPos)
    {
        final GeneAnnotation geneAnnotation = trans.parent();

        long position = geneAnnotation.position();

        if(requireHigherBreakendPos)
        {
            // factor in any uncertainty around the precise breakend, eg from homology
            if(geneAnnotation.variant() != null)
                position += geneAnnotation.isStart() ? geneAnnotation.variant().start().endOffset() : geneAnnotation.variant().end().endOffset();

            return (position >= rnaPosition);
        }
        else
        {
            if(geneAnnotation.variant() != null)
                position += geneAnnotation.isStart() ? geneAnnotation.variant().start().startOffset() : geneAnnotation.variant().end().startOffset();

            return (position <= rnaPosition);
        }
    }

    public void setRnaFusionData(final RnaFusionData rnaFusion)
    {
        int[] transUpExonData = mGeneTranscriptCollection.getExonData(rnaFusion.GeneUp, rnaFusion.PositionUp);
        rnaFusion.setExonUpRank(transUpExonData[EXON_RANK_MIN], transUpExonData[EXON_RANK_MAX]);

        transUpExonData = mGeneTranscriptCollection.getExonData(rnaFusion.GeneDown, rnaFusion.PositionDown);
        rnaFusion.setExonDownRank(transUpExonData[EXON_RANK_MIN], transUpExonData[EXON_RANK_MAX]);
    }

    public void onCompleted()
    {
        closeBufferedWriter(mFusionWriter);
        closeBufferedWriter(mRnaWriter);
    }


}
