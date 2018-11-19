package com.hartwig.hmftools.svanalysis.analysis;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svannotation.analysis.SvFusionAnalyser.FUSION_PAIRS_CSV;
import static com.hartwig.hmftools.svannotation.analysis.SvFusionAnalyser.PROMISCUOUS_FIVE_CSV;
import static com.hartwig.hmftools.svannotation.analysis.SvFusionAnalyser.PROMISCUOUS_THREE_CSV;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.SvGeneTranscriptCollection;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvVarData;
import com.hartwig.hmftools.svannotation.analysis.SvDisruptionAnalyser;
import com.hartwig.hmftools.svannotation.analysis.SvFusionAnalyser;

import org.apache.commons.cli.CommandLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class FusionDisruptionAnalyser
{
    private SvFusionAnalyser mFusionFinder;
    private SvDisruptionAnalyser mDisruptionFinder;

    private String mSampleId;
    private String mOutputDir;
    private boolean mUseCombinedOutput;
    private SvGeneTranscriptCollection mSvGeneTranscriptCollection;

    private List<GeneFusion> mGeneFusions;
    private List<GeneDisruption> mGeneDisruptions;

    private BufferedWriter mFusionWriter;

    private static final Logger LOGGER = LogManager.getLogger(FusionDisruptionAnalyser.class);

    public FusionDisruptionAnalyser()
    {
        mFusionFinder = null;
        mDisruptionFinder = null;
        mSvGeneTranscriptCollection = new SvGeneTranscriptCollection();

        mGeneFusions = Lists.newArrayList();
        mGeneDisruptions = Lists.newArrayList();
        mFusionWriter= null;
        mOutputDir = "";
        mUseCombinedOutput = false;
    }

    public boolean loadFusionReferenceData(final CommandLine cmdLineArgs, final String outputDir, boolean useCombinedOutput)
    {
        try
        {
            KnownFusionsModel knownFusionsModel = KnownFusionsModel.fromInputStreams(
                    new FileInputStream(cmdLineArgs.getOptionValue(FUSION_PAIRS_CSV)),
                    new FileInputStream(cmdLineArgs.getOptionValue(PROMISCUOUS_FIVE_CSV)),
                    new FileInputStream(cmdLineArgs.getOptionValue(PROMISCUOUS_THREE_CSV)));

            mFusionFinder = new SvFusionAnalyser(knownFusionsModel);

        }
        catch(IOException e)
        {
            LOGGER.error("failed to load known fusion files");
            return false;
        }

        mOutputDir = outputDir;
        mUseCombinedOutput = useCombinedOutput;

        return true;
    }

    public void loadSvGeneTranscriptData(final String sampleId, final String sampleDataPath)
    {
        mSampleId = sampleId;
        mSvGeneTranscriptCollection.setDataPath(sampleDataPath);
        mSvGeneTranscriptCollection.loadSampleGeneTranscripts(sampleId);
    }

    private final List<GeneAnnotation> getSvGenesList(final SvVarData var, boolean isStart)
    {
        final List<GeneAnnotation> genesList = mSvGeneTranscriptCollection.getSvIdGeneTranscriptsMap().get(var.dbId());

        if(genesList == null)
            return Lists.newArrayList();

        if(isStart)
            return genesList.stream().filter(GeneAnnotation::isStart).collect(Collectors.toList());
        else
            return genesList.stream().filter(GeneAnnotation::isEnd).collect(Collectors.toList());
    }

    // private static String CHECK_VAR_ID = "";
    private static String CHECK_VAR_ID = "468339";

    public void findFusions(final List<SvVarData> allVariants, final List<SvCluster> clusters, final List<SvVarData> svList)
    {
        if(mSampleId.isEmpty() || mFusionFinder == null)
            return;

        final Map<Integer, List<GeneAnnotation>> svGenesMap = mSvGeneTranscriptCollection.getSvIdGeneTranscriptsMap();

        if(svGenesMap.isEmpty())
            return;

        List<StructuralVariantData> svDataList = allVariants.stream().map(SvVarData::getSvData).collect(Collectors.toList());

        mSvGeneTranscriptCollection.setSvData(svDataList);

        // always report SVs by themselves
        for(final SvVarData var : svList)
        {
            if(var.isNullBreakend())
                continue;

            if(var.id().equals(CHECK_VAR_ID))
            {
                LOGGER.debug("specific var({})", var.posId());
            }

            List<GeneAnnotation> breakendGenes1 = getSvGenesList(var, true);
            List<GeneAnnotation> breakendGenes2 = getSvGenesList(var, false);

            checkFusions(breakendGenes1, breakendGenes2, var.getCluster());
        }

        boolean checkClusters = true;

        if(checkClusters)
        {
            // for now only consider simple SVs and resolved small clusters
            for (final SvCluster cluster : clusters)
            {
                List<GeneAnnotation> breakendGenes1 = Lists.newArrayList();
                List<GeneAnnotation> breakendGenes2 = Lists.newArrayList();

                /*
                if(cluster.getId() == 651)
                {
                    LOGGER.debug("specific cluster");
                }
                */

                if (cluster.getCount() == 1)
                {
                    // checked above already
                    continue;
                }
                else if (cluster.isResolved() && cluster.isFullyChained() && cluster.isConsistent() && cluster.getTypeCount(SGL) == 0)
                {
                    final SvChain completeChain = cluster.getChains().get(0);
                    final SvVarData startVar = completeChain.getFirstSV();
                    final SvVarData endVar = completeChain.getLastSV();

                    breakendGenes1 = getSvGenesList(startVar, completeChain.firstLinkOpenOnStart());
                    breakendGenes2 = getSvGenesList(endVar, completeChain.lastLinkOpenOnStart());
                }

                checkFusions(breakendGenes1, breakendGenes2, cluster);
            }
        }
    }

    private void checkFusions(List<GeneAnnotation> breakendGenes1, List<GeneAnnotation> breakendGenes2, final SvCluster cluster)
    {
        if(breakendGenes1.isEmpty() || breakendGenes2.isEmpty())
            return;

        List<GeneFusion> fusions = mFusionFinder.findFusions(breakendGenes1, breakendGenes2);

        if(fusions.isEmpty())
            return;

        // for now only log reportable fusions
        fusions = fusions.stream().filter(GeneFusion::reportable).collect(Collectors.toList());

        if(LOGGER.isDebugEnabled())
        {
            for (final GeneFusion fusion : fusions)
            {
                if (fusion.reportable())
                {
                    final Transcript upstream = fusion.upstreamLinkedAnnotation();
                    final Transcript downstream = fusion.downstreamLinkedAnnotation();
                    final GeneAnnotation upGene = upstream.parent();
                    final GeneAnnotation downGene = downstream.parent();

                    LOGGER.debug("sample({}) fusion: up({} {} {} {}) upSV({}: {}:{}:{} start={} strand={}) up({} {} {} {}) upSV({}: {}:{}:{} start={} strand={})",
                            mSampleId, upstream.geneName(), upstream.transcriptId(), upstream.regionType(), upstream.codingType(),
                            upGene.id(), upGene.chromosome(), upGene.position(), upGene.orientation(), upGene.isStart(), upGene.strand(),
                            downstream.geneName(), downstream.transcriptId(), downstream.regionType(), downstream.codingType(),
                            downGene.id(), downGene.chromosome(), downGene.position(), downGene.orientation(), downGene.isStart(), downGene.strand());
                }
            }
        }

        writeFusions(fusions, cluster);
    }

    private void writeFusions(final List<GeneFusion> fusions, final SvCluster cluster)
    {
        if(fusions.isEmpty())
            return;

        try
        {
            BufferedWriter writer = null;

            if(mFusionWriter == null)
            {
                String outputFilename = mOutputDir;

                if (!outputFilename.endsWith("/"))
                    outputFilename += "/";

                if(mUseCombinedOutput)
                    outputFilename += "FUSIONS.csv";
                else
                    outputFilename += mSampleId + "_" + "sv_fusions.csv";

                Path outputFile = Paths.get(outputFilename);

                mFusionWriter = Files.newBufferedWriter(outputFile, StandardOpenOption.CREATE);
                writer = mFusionWriter;

                writer.write("SampleId,Reportable,PrimarySource,ClusterId,ClusterCount,ResolvedType,SameVar");
                writer.write(",StartSvId,StartChr,StartPos,StartOrient,StartType");
                writer.write(",StartGene,StartTranscript,StartStrand,StartRegionType,StartCodingType,StartExon,StartPhase,StartCodingBases,StartTotalCodingBases,StartCodingStart,StartCodingEnd");
                writer.write(",EndSvId,EndChr,EndPos,EndOrient,EndType");
                writer.write(",EndGene,EndTranscript,EndStrand,EndRegionType,EndCodingType,EndExon,EndPhase,EndCodingBases,EndTotalCodingBases,EndCodingStart,EndCodingEnd");
                writer.newLine();
            }
            else
            {
                writer = mFusionWriter;
            }

            for(final GeneFusion fusion : fusions)
            {
                final Transcript startTrans = fusion.upstreamLinkedAnnotation();
                final Transcript endTrans = fusion.downstreamLinkedAnnotation();

                final GeneAnnotation startVar = startTrans.parent();
                final GeneAnnotation endVar = endTrans.parent();

                boolean sameVar = startVar.id() == endVar.id();

                writer.write(String.format("%s,%s,%s,%d,%d,%s,%s",
                        mSampleId, fusion.reportable(), fusion.primarySource(),
                        cluster.getId(), cluster.getUniqueSvCount(), cluster.getResolvedType(), sameVar));

                // write upstream SV, transcript and exon info
                writer.write(
                        String.format(",%d,%s,%d,%d,%s",
                                startVar.id(), startVar.chromosome(), startVar.position(), startVar.orientation(), startVar.type()));

                writer.write(
                        String.format(",%s,%s,%d,%s,%s,%d,%d,%d,%d,%d,%d",
                                startTrans.parent().geneName(), startTrans.transcriptId(),
                                startTrans.parent().strand(), startTrans.regionType(), startTrans.codingType(),
                                startTrans.exonUpstream(), startTrans.exonUpstreamPhase(),
                                startTrans.codingBases(), startTrans.totalCodingBases(),
                                startTrans.codingStart() != null ? startTrans.codingStart() : 0,
                                startTrans.codingEnd() != null ? startTrans.codingEnd() : 0));

                writer.write(
                        String.format(",%d,%s,%d,%d,%s",
                                endVar.id(), endVar.chromosome(), endVar.position(), endVar.orientation(), endVar.type()));

                writer.write(
                        String.format(",%s,%s,%d,%s,%s,%d,%d,%d,%d,%d,%d",
                                endTrans.parent().geneName(), endTrans.transcriptId(),
                                endTrans.parent().strand(), endTrans.regionType(), endTrans.codingType(),
                                endTrans.exonUpstream(), endTrans.exonUpstreamPhase(),
                                endTrans.codingBases(), endTrans.totalCodingBases(),
                                endTrans.codingStart() != null ? endTrans.codingStart() : 0,
                                endTrans.codingEnd() != null ? endTrans.codingEnd() : 0));

                writer.newLine();
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing fusions: {}", e.toString());
        }
    }

    public void close()
    {
        try
        {
            if(mFusionWriter != null)
                mFusionWriter.close();
        }
        catch (IOException e)
        {

        }
    }
}
