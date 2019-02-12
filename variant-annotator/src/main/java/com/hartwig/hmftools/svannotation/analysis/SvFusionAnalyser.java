package com.hartwig.hmftools.svannotation.analysis;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation.isDownstream;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation.isUpstream;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_3P_PROM;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_5P_PROM;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_BOTH_PROM;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_KNOWN;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_NONE;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.EXON_RANK_MAX;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.EXON_RANK_MIN;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
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

    public static final String SAMPLE_RNA_FILE = "sample_rna_file";

    private static final int EXON_THRESHOLD = 1;

    private KnownFusionsModel mKnownFusionsModel;

    private Map<String, List<RnaFusionData>> mSampleRnaData;
    private SvGeneTranscriptCollection mGeneTranscriptCollection;
    private List<String> mProteinsRequiredKept;
    private List<String> mProteinsRequiredLost;

    private BufferedWriter mFusionWriter;
    private BufferedWriter mRnaWriter;

    private static final Logger LOGGER = LogManager.getLogger(SvFusionAnalyser.class);

    public SvFusionAnalyser(final CommandLine cmd, final SvGeneTranscriptCollection geneTranscriptCollection)
    {
        mGeneTranscriptCollection = geneTranscriptCollection;

        mKnownFusionsModel = null;

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

        mFusionWriter = null;
        mRnaWriter = null;
        mSampleRnaData = Maps.newHashMap();

        if (cmd.hasOption(SAMPLE_RNA_FILE))
        {
            loadSampleRnaData(cmd.getOptionValue(SAMPLE_RNA_FILE));
        }

        mProteinsRequiredKept = Lists.newArrayList();
        mProteinsRequiredLost = Lists.newArrayList();
        setRequiredProteins();
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(FUSION_PAIRS_CSV, true, "Path towards a CSV containing white-listed gene fusion pairs.");
        options.addOption(PROMISCUOUS_FIVE_CSV, true, "Path towards a CSV containing white-listed promiscuous 5' genes.");
        options.addOption(PROMISCUOUS_THREE_CSV, true, "Path towards a CSV containing white-listed promiscuous 3' genes.");
        options.addOption(SAMPLE_RNA_FILE, true, "Sample RNA data to match");
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

    private static int SPECIFIC_VAR_ID = -1;
    // private static int SPECIFIC_VAR_ID = 13699019;

    // private static String SPEC_GENE_1 = "ZMYND12";
    // private static String SPEC_GENE_2 = "ABL1";
    private static String SPEC_GENE_1 = "";
    private static String SPEC_GENE_2 = "";

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

                if((startGene.GeneName.equals(SPEC_GENE_1) && endGene.GeneName.equals(SPEC_GENE_2))
                || (startGene.GeneName.equals(SPEC_GENE_2) && endGene.GeneName.equals(SPEC_GENE_1)))
                {
                    LOGGER.debug("gene match: {} and {}", startGene.GeneName, endGene.GeneName);
                }

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

                        if(downstreamTrans.exonMax() == 1)
                            continue;

                        if(!upstreamTrans.isDisruptive())
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
                                addFusion(potentialFusions, upstreamTrans, downstreamTrans);
                            }
                        }
                        else
                        {
                            // just check for a phasing match
                            if (upstreamTrans.exonUpstreamPhase() == downstreamTrans.exonDownstreamPhase())
                            {
                                addFusion(potentialFusions, upstreamTrans, downstreamTrans);
                            }
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

    private void addFusion(List<GeneFusion> fusions, final Transcript upstreamTrans, final Transcript downstreamTrans)
    {
        //LOGGER.debug("adding fusion between start SV({}) trans({}) and end SV({}) trans({})",
        //        startTrans.parent().id(), startTrans.toString(), endTrans.parent().id(), endTrans.toString());

        fusions.add(new GeneFusion(
                upstreamTrans,
                downstreamTrans,
                mKnownFusionsModel.primarySource(upstreamTrans.parent().synonyms(), downstreamTrans.parent().synonyms()),
                false));
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
        if(t1.isIntronic() && t2.isIntronic() && t1.exonUpstream() == t2.exonUpstream())
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
        && reportableFusion.upstreamTrans().exonUpstreamPhase() == -1)
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

            long length = downTrans.isCoding() ? downTrans.calcCodingBases(false) : downTrans.exonMax();

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

            length = upTrans.isCoding() ? upTrans.calcCodingBases(true) : upTrans.exonMax();
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
                && downTrans.exonDownstream() - upTrans.exonUpstream() > EXON_THRESHOLD);

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

    public void writeFusions(final List<GeneFusion> fusions, final String outputDir, final String sampleId,  final String clusterInfo, boolean hasMultipleSamples)
    {
        try
        {
            if(mFusionWriter == null)
            {
                String outputFilename = outputDir;

                if (!outputFilename.endsWith(File.separator))
                    outputFilename += File.separator;

                if(hasMultipleSamples)
                    outputFilename += "FUSIONS.csv";
                else
                    outputFilename += sampleId + "_fusions.csv";

                mFusionWriter = createBufferedWriter(outputFilename, false);

                mFusionWriter.write("SampleId,Reportable,KnownType,PrimarySource,ClusterId,ClusterCount,ResolvedType");

                mFusionWriter.write(",SvIdUp,ChrUp,PosUp,OrientUp,TypeUp,PloidyUp,GeneUp,ChrBandUp,TranscriptUp,StrandUp,RegionTypeUp,CodingTypeUp");
                mFusionWriter.write(",ExonUp,PhaseUp,ExonMaxUp,DisruptiveUp,ExactBaseUp,CodingBasesUp,TotalCodingUp");
                mFusionWriter.write(",CodingStartUp,CodingEndUp,TransStartUp,TransEndUp,DistancePrevUp,CanonicalUp,BiotypeUp");

                mFusionWriter.write(",SvIdDown,ChrDown,PosDown,OrientDown,TypeDown,PloidyDown,GeneDown,ChrBandDown,TranscriptDown,StrandDown,RegionTypeDown,CodingTypeDown");
                mFusionWriter.write(",ExonDown,PhaseDown,ExonMaxDown,DisruptiveDown,ExactBaseDown,CodingBasesDown,TotalCodingDown");
                mFusionWriter.write(",CodingStartDown,CodingEndDown,TransStartDown,TransEndDown,DistancePrevDown,CanonicalDown,BiotypeDown");

                mFusionWriter.write(",ProteinsKept,ProteinsLost");
                mFusionWriter.newLine();
            }

            BufferedWriter writer = mFusionWriter;

            for(final GeneFusion fusion : fusions)
            {
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
                                upTrans.exonUpstream(), upTrans.exonUpstreamPhase(), upTrans.exonMax(), upTrans.isDisruptive()));
                writer.write(
                        String.format(",%d,%d,%d,%d,%d,%d,%d,%d,%s,%s",
                                upTrans.exactCodingBase(), upTrans.calcCodingBases(true), upTrans.totalCodingBases(),
                                upTrans.codingStart(), upTrans.codingEnd(), upTrans.transcriptStart(), upTrans.transcriptEnd(),
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
                                downTrans.exonDownstream(), downTrans.exonDownstreamPhase(), downTrans.exonMax(), downTrans.isDisruptive()));

                writer.write(
                        String.format(",%d,%d,%d,%d,%d,%d,%d,%d,%s,%s,%s,%s",
                                downTrans.exactCodingBase(), downTrans.calcCodingBases(false), downTrans.totalCodingBases(),
                                downTrans.codingStart(), downTrans.codingEnd(), downTrans.transcriptStart(), downTrans.transcriptEnd(),
                                downTrans.exonDistanceUp(), downTrans.isCanonical(), downTrans.bioType(),
                                downTrans.getProteinFeaturesKept(), downTrans.getProteinFeaturesLost()));

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
                                    transUp.parent().orientation(), transUp.parent().Strand, transUp.parent().type()));

                    writer.write(
                            String.format(",%s,%s,%s,%d,%d,%s,%d,%d",
                                    transUp.StableId, transUp.regionType(), transUp.codingType(),
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
                                    transDown.parent().orientation(), transDown.parent().Strand, transDown.parent().type()));

                    writer.write(
                            String.format(",%s,%s,%s,%d,%d,%s,%d,%d",
                                    transDown.StableId, transDown.regionType(), transDown.codingType(),
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

    private static String SPECIFIC_FUSION = "";
    // private static String SPECIFIC_FUSION = "TRPM7--RP11-184D12.1";

    private void matchRnaFusion(final RnaFusionData rnaFusion, final List<GeneFusion> fusions, final List<StructuralVariantAnnotation> annotations)
    {
        // find the min and max exon rank for every transcript matching the RNA positions
        setRnaFusionData(rnaFusion);

        if(rnaFusion.Name.equals(SPECIFIC_FUSION))
        {
            LOGGER.debug("specific RNA fusion: {}", rnaFusion.Name);
        }

        // first check for a fusion found
        for(final GeneFusion fusion : fusions)
        {
            final Transcript upTrans = fusion.upstreamTrans();
            final Transcript downTrans = fusion.downstreamTrans();

            if(!upTrans.geneName().equals(rnaFusion.GeneUp)  || !downTrans.geneName().equals(rnaFusion.GeneDown))
                continue;

            // breakends must be downstream of the upstream RNA exons, and upstream of the downstream RNA exon
            if(!withinPositionRange(upTrans, rnaFusion.PositionUp, upTrans.parent().Strand == 1))
                continue;

            if(!withinPositionRange(downTrans, rnaFusion.PositionDown, downTrans.parent().Strand != 1))
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
                Transcript svTransUp = null;
                Transcript svTransDown = null;

                for (final GeneAnnotation breakendGene : svAnnotation.annotations())
                {
                    for(final Transcript transcript : breakendGene.transcripts())
                    {
                        if(isUpstream(breakendGene) && transcript.geneName().equals(rnaFusion.GeneUp)
                        && withinPositionRange(transcript, rnaFusion.PositionUp, transcript.parent().Strand == 1))
                        {
                            transUp = transcript;
                            svTransUp = transcript;
                        }

                        if(isDownstream(breakendGene) && transcript.geneName().equals(rnaFusion.GeneDown)
                        && withinPositionRange(transcript, rnaFusion.PositionDown, transcript.parent().Strand != 1))
                        {
                            transDown= transcript;
                            svTransDown = transcript;
                        }

                        if(svTransDown != null && svTransUp != null)
                            break;
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
