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
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;
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
            mFusionFinder.setIncludePossibles(true);

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

    private void setSvGenesList(final SvVarData var)
    {
        final List<GeneAnnotation> genesList = mSvGeneTranscriptCollection.getSvIdGeneTranscriptsMap().get(var.dbId());

        List<GeneAnnotation> startGenes = var.getGenesList(true);
        List<GeneAnnotation> endGenes = var.getGenesList(false);

        for(GeneAnnotation gene : genesList)
        {
            gene.setSvData(var.getSvData());

            if(gene.isStart())
                startGenes.add(gene);
            else
                endGenes.add(gene);
        }
    }

    private static String CHECK_VAR_ID = "";
    // private static String CHECK_VAR_ID = "527632";

    public void findFusions(final List<SvVarData> svList, final List<SvCluster> clusters, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        if(mSampleId.isEmpty() || mFusionFinder == null)
            return;

        final Map<Integer, List<GeneAnnotation>> svGenesMap = mSvGeneTranscriptCollection.getSvIdGeneTranscriptsMap();

        if(svGenesMap.isEmpty())
            return;

        List<GeneFusion> allFusions = Lists.newArrayList();

        // always report SVs by themselves
        for(final SvVarData var : svList)
        {
            if(var.isReplicatedSv())
                continue;

            // cache transcript info against each SV
            setSvGenesList(var);

            if(var.isNullBreakend())
                continue;

            if(var.id().equals(CHECK_VAR_ID))
            {
                LOGGER.debug("specific var({})", var.posId());
            }

            checkFusions(var.getGenesList(true), var.getGenesList(false), var.getCluster());
        }

        boolean checkClusters = false;
        int maxClusterSize = 10;

        if(checkClusters)
        {
            // for now only consider simple SVs and resolved small clusters
            for (final SvCluster cluster : clusters)
            {
                /*
                if(cluster.getId() == 651)
                {
                    LOGGER.debug("specific cluster");
                }
                */

                if (cluster.getCount() == 1) // simple clusters already checked
                    continue;

                if(cluster.hasReplicatedSVs() || !cluster.isFullyChained() || cluster.getTypeCount(SGL) > 0)
                    continue;

                if(cluster.getUniqueSvCount() > maxClusterSize)
                    continue;

                // test every breakend with every other one in the chain since the fusion could be formed at any point
                GeneFusion fusion = findChainedFusion(cluster.getChains().get(0));
            }
        }
    }

    private GeneFusion findChainedFusion(final SvChain chain)
    {
        for(int index1 = 0; index1 < chain.getSvList().size()-1; ++index1)
        {
            final SvVarData var1 = chain.getSvList().get(index1);

            for(int index2 = index1+1; index2 < chain.getSvList().size(); ++index2)
            {
                final SvVarData var2 = chain.getSvList().get(index2);


            }

        }

        /*
        for(final SvLinkedPair linkedPair : chain.getLinkedPairs())
        {
            List<GeneAnnotation> breakendGenes1 = null;
            List<GeneAnnotation> breakendGenes2 = null;
        }


            checkFusions(breakendGenes1, breakendGenes2, cluster);

            final SvVarData startVar = completeChain.getFirstSV();
            final SvVarData endVar = completeChain.getLastSV();

            breakendGenes1 = getSvGenesList(startVar, completeChain.firstLinkOpenOnStart());
            breakendGenes2 = getSvGenesList(endVar, completeChain.lastLinkOpenOnStart());


        }
        */

        return null;
    }

    private boolean isFusionDisrupted(final GeneFusion fusion, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        return false;
    }

    private void checkFusions(List<GeneAnnotation> breakendGenes1, List<GeneAnnotation> breakendGenes2, final SvCluster cluster)
    {
        if (breakendGenes1.isEmpty() || breakendGenes2.isEmpty())
            return;

        List<GeneFusion> fusions = mFusionFinder.findFusions(breakendGenes1, breakendGenes2);

        if (fusions.isEmpty())
            return;

        // fusions = fusions.stream().filter(GeneFusion::reportable).collect(Collectors.toList()); // restrict to reportable fusions

        if(LOGGER.isDebugEnabled())
        {
            for (final GeneFusion fusion : fusions)
            {
                if (fusion.reportable())
                {
                    final Transcript upstream = fusion.upstreamTrans();
                    final Transcript downstream = fusion.downstreamTrans();
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

                writer.write("SampleId,Reportable,PrimarySource,ClusterId,ClusterCount,ResolvedType,PhaseMatched");
                writer.write(",SvIdUp,ChrUp,PosUp,OrientUp,TypeStart,GeneUp,TranscriptUp,StrandUp,RegionTypeUp,CodingTypeUp");
                writer.write(",ExonUp,PhaseUp,ExactBaseUp,CodingBasesUp,TotalCodingUp,ExonMaxUp,CodingStartUp,CodingEndUp");
                writer.write(",SvIdDown,ChrDown,PosDown,OrientDown,TypeDown,GeneDown,TranscriptDown,StrandDown,RegionTypeDown,CodingTypeDown");
                writer.write(",ExonDown,PhaseDown,ExactBaseDown,CodingBasesDown,TotalCodingDown,ExonMaxDown,CodingStartDown,CodingEndDown");
                writer.newLine();
            }
            else
            {
                writer = mFusionWriter;
            }

            for(final GeneFusion fusion : fusions)
            {
                final Transcript startTrans = fusion.upstreamTrans();
                final Transcript endTrans = fusion.downstreamTrans();

                final GeneAnnotation startVar = startTrans.parent();
                final GeneAnnotation endVar = endTrans.parent();

                writer.write(String.format("%s,%s,%s,%d,%d,%s,%s",
                        mSampleId, fusion.reportable(), fusion.primarySource(),
                        cluster.getId(), cluster.getUniqueSvCount(), cluster.getResolvedType(), fusion.isPhaseMatch()));

                // write upstream SV, transcript and exon info
                writer.write(
                        String.format(",%d,%s,%d,%d,%s",
                                startVar.id(), startVar.chromosome(), startVar.position(), startVar.orientation(), startVar.type()));

                writer.write(
                        String.format(",%s,%s,%d,%s,%s",
                                startTrans.parent().geneName(), startTrans.transcriptId(),
                                startTrans.parent().strand(), startTrans.regionType(), startTrans.codingType()));

                writer.write(
                        String.format(",%d,%d,%d,%d,%d,%d,%d,%d",
                                startTrans.exonUpstream(), startTrans.exonUpstreamPhase(), startTrans.exactCodingBase(),
                                startTrans.codingBases(), startTrans.totalCodingBases(), startTrans.exonMax(),
                                startTrans.codingStart() != null ? startTrans.codingStart() : 0,
                                startTrans.codingEnd() != null ? startTrans.codingEnd() : 0));

                writer.write(
                        String.format(",%d,%s,%d,%d,%s",
                                endVar.id(), endVar.chromosome(), endVar.position(), endVar.orientation(), endVar.type()));

                writer.write(
                        String.format(",%s,%s,%d,%s,%s",
                                endTrans.parent().geneName(), endTrans.transcriptId(),
                                endTrans.parent().strand(), endTrans.regionType(), endTrans.codingType()));

                writer.write(
                        String.format(",%d,%d,%d,%d,%d,%d,%d,%d",
                                endTrans.exonDownstream(), endTrans.exonDownstreamPhase(), endTrans.exactCodingBase(),
                                endTrans.codingBases(), endTrans.totalCodingBases(), endTrans.exonMax(),
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
