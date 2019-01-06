package com.hartwig.hmftools.svannotation.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
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
import com.hartwig.hmftools.svannotation.EnsemblDAO;
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
    // private static int SPECIFIC_VAR_ID = 4740717;

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
                .filter(GeneFusion::isPhaseMatch)
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

    private static boolean isUpstream(GeneAnnotation gene)
    {
        return gene.strand() * gene.orientation() > 0;
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

    public void onCompleted()
    {
        closeBufferedWriter(mFusionWriter);
        closeBufferedWriter(mRnaWriter);
    }

    public final Map<String, List<RnaFusionData>> getSampleRnaData() { return mSampleRnaData; }
    public final List<RnaFusionData> getSampleRnaData(final String sampleId) { return mSampleRnaData.get(sampleId); }

    private static int COL_SAMPLEID = 0;
    private static int COL_NAME = 1;
    private static int COL_GENE_UP = 5;
    private static int COL_POS_UP = 8;
    private static int COL_ORIENT_UP = 9;
    private static int COL_GENE_DOWN = 10;
    private static int COL_POS_DOWN = 13;
    private static int COL_ORIENT_DOWN = 14;

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
                        items[COL_GENE_UP], items[COL_GENE_DOWN],
                        Long.parseLong(items[COL_POS_UP]), Long.parseLong(items[COL_POS_DOWN]),
                        Byte.parseByte(items[COL_ORIENT_UP]), Byte.parseByte(items[COL_ORIENT_DOWN]));

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

                mRnaWriter.write("SampleId,FusionName,GeneUp,GeneDown");
                mRnaWriter.write(",SvFusionMatch,Reportable,SvIdUp,SvIdDown,ExonRankUp,ExonRankDown");
                mRnaWriter.write(",TransExonMinUp,TransExonMaxUp,TransExonMinDown,TransExonMaxDown");
                mRnaWriter.newLine();
            }

            BufferedWriter writer = mRnaWriter;

            for(final RnaFusionData rnaFusion : rnaFusionList)
            {
                // first check for a fusion found
                GeneFusion fusionMatch = null;

                for(final GeneFusion fusion : fusions)
                {
                    if(fusion.upstreamTrans().geneName().equals(rnaFusion.GeneUp) && fusion.downstreamTrans().geneName().equals(rnaFusion.GeneDown))
                    {
                        fusionMatch = fusion;
                        break;
                    }
                }

                GeneAnnotation upGeneAnnotation = null;
                GeneAnnotation downGeneAnnotation = null;

                if(fusionMatch != null)
                {
                    upGeneAnnotation = fusionMatch.upstreamTrans().parent();
                    downGeneAnnotation = fusionMatch.downstreamTrans().parent();
                }
                else
                {
                    // check breakend annotations for any match on the gene
                    for(final StructuralVariantAnnotation svAnnotation : svAnnotations)
                    {
                        List<GeneAnnotation> breakendGenes = Lists.newArrayList(svAnnotation.start());
                        breakendGenes.addAll(svAnnotation.end());

                        for (final GeneAnnotation breakendGene : breakendGenes)
                        {
                            if (isUpstream(breakendGene) && breakendGene.geneName().equals(rnaFusion.GeneUp))
                            {
                                upGeneAnnotation = breakendGene;
                            }
                            else if (!isUpstream(breakendGene) && breakendGene.geneName().equals(rnaFusion.GeneDown))
                            {
                                downGeneAnnotation = breakendGene;
                            }

                            if(upGeneAnnotation != null && downGeneAnnotation != null)
                                break;
                        }
                    }
                }

                // find the min and max exon rank for every transcript matching the RNA positions
                setRnaFusionData(rnaFusion);

                writer.write(
                        String.format(",%s,%s,%s,%s",
                                sampleId, rnaFusion.Name, rnaFusion.GeneUp, rnaFusion.GeneDown));

                if(fusionMatch != null)
                {
                    writer.write(String.format(",%s,%s,%d,%d,%d,%d",
                            "true", fusionMatch.reportable(),
                            fusionMatch.upstreamTrans().parent().id(), fusionMatch.downstreamTrans().parent().id(),
                            fusionMatch.upstreamTrans().exonUpstream(), fusionMatch.downstreamTrans().exonDownstream()));
                }
                else
                {
                    writer.write(String.format(",false,false,%d,%d,-1,-1",
                            upGeneAnnotation != null ? upGeneAnnotation.id() : -1,
                            downGeneAnnotation != null ? downGeneAnnotation.id() : -1));
                }

                writer.write(String.format(",%d,%d,%d,%d",
                        rnaFusion.exonMinRankUp(), rnaFusion.exonMaxRankUp(), rnaFusion.exonMinRankDown(), rnaFusion.exonMaxRankDown()));

                writer.newLine();

            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing RNA match data: {}", e.toString());
        }

    }

    private static int EXON_RANK_MIN = 0;
    private static int EXON_RANK_MAX = 1;

    public void setRnaFusionData(final RnaFusionData rnaFusion)
    {
        int[] transUpExonData = getExonData(rnaFusion.GeneUp, rnaFusion.PositionUp);
        rnaFusion.setExonUpRank(transUpExonData[EXON_RANK_MIN], transUpExonData[EXON_RANK_MAX]);

        transUpExonData = getExonData(rnaFusion.GeneDown, rnaFusion.PositionDown);
        rnaFusion.setExonDownRank(transUpExonData[EXON_RANK_MIN], transUpExonData[EXON_RANK_MAX]);
    }

    private int[] getExonData(final String geneName, long exonPosition)
    {
        int[] exonData = new int[EXON_RANK_MAX+1];

        final EnsemblGeneData geneData = mGeneTranscriptCollection.getGeneData(geneName);

        if(geneData == null)
            return exonData;

        final List<TranscriptExonData> transExonDataList = mGeneTranscriptCollection.getTransExonData(geneData.GeneId);

        if(transExonDataList == null)
            return exonData;

        int minRank = -1;
        int maxRank = -1;
        for(final TranscriptExonData transExonData : transExonDataList)
        {
            int exonRank = transExonData.ExonRank;
            // final String transcriptStableId = transcriptData.get(TRANSCRIPT.STABLE_ID);
            // final UInteger transcriptId = transcriptData.get(TRANSCRIPT.TRANSCRIPT_ID);

            minRank = minRank == -1 ? exonRank : min(exonRank, minRank);
            maxRank = max(exonRank, maxRank);
        }

        exonData[EXON_RANK_MIN] = minRank;
        exonData[EXON_RANK_MAX] = maxRank;
        return exonData;

    }

}
