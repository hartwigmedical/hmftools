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
    private SvGeneTranscriptCollection mSvGeneTranscriptCollection;

    private List<GeneFusion> mGeneFusions;
    private List<GeneDisruption> mGeneDisruptions;

    private static final Logger LOGGER = LogManager.getLogger(FusionDisruptionAnalyser.class);

    public FusionDisruptionAnalyser()
    {
        mFusionFinder = null;
        mDisruptionFinder = null;
        mSvGeneTranscriptCollection = new SvGeneTranscriptCollection();

        mGeneFusions = Lists.newArrayList();
        mGeneDisruptions = Lists.newArrayList();
    }

    public boolean loadFusionReferenceData(final CommandLine cmdLineArgs)
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

    public void findFusions(final List<SvVarData> allVariants, final List<SvCluster> clusters, final String outputDir)
    {
        if(mSampleId.isEmpty() || mFusionFinder == null)
            return;

        final Map<Integer, List<GeneAnnotation>> svGenesMap = mSvGeneTranscriptCollection.getSvIdGeneTranscriptsMap();

        if(svGenesMap.isEmpty())
            return;

        List<StructuralVariantData> svDataList = allVariants.stream().map(SvVarData::getSvData).collect(Collectors.toList());

        mSvGeneTranscriptCollection.setSvData(svDataList);

        List<GeneAnnotation> breakendGenes1 = Lists.newArrayList();
        List<GeneAnnotation> breakendGenes2 = Lists.newArrayList();

        // for now only consider simple SVs and resolved small clusters
        for(final SvCluster cluster : clusters)
        {
            breakendGenes1.clear();
            breakendGenes2.clear();

            if(cluster.getCount() == 1)
            {
                final SvVarData var = cluster.getSVs().get(0);

                breakendGenes1 = getSvGenesList(var, true);
                breakendGenes2 = getSvGenesList(var, false);
            }
            else if(cluster.isResolved() && cluster.isFullyChained() && cluster.isConsistent() && cluster.getTypeCount(SGL) == 0)
            {
                final SvChain completeChain = cluster.getChains().get(0);
                final SvVarData startVar = completeChain.getFirstSV();
                final SvVarData endVar = completeChain.getLastSV();

                breakendGenes1 = getSvGenesList(startVar, completeChain.firstLinkOpenOnStart());
                breakendGenes2 = getSvGenesList(endVar, completeChain.lastLinkOpenOnStart());
            }

            if(breakendGenes1.isEmpty() || breakendGenes2.isEmpty())
                continue;

            List<GeneFusion> fusions = mFusionFinder.findFusions(breakendGenes1, breakendGenes2);

        }

        // writeFusions
    }

    private void writeFusions(final List<GeneFusion> fusions, final String outputDir)
    {
        if(fusions.isEmpty())
            return;

        String outputFilename = outputDir;

        if(!outputFilename.endsWith("/"))
            outputFilename += "/";

        outputFilename += mSampleId + "_" + "sv_fusions.csv";

        try
        {
            Path outputFile = Paths.get(outputFilename);

            BufferedWriter writer = Files.newBufferedWriter(outputFile, StandardOpenOption.CREATE);

            writer.write("SampleId");
            writer.write(",StartSvId,StartChr,StartPos,StartOrient,StartType");
            writer.write(",StartGene,StartTranscript,StartRegionType,StartExon,StartPhase,StartCodingBases,StartTotalCodingBases,StartCodingStart,StartCodingEnd");
            writer.write(",EndSvId,EndChr,EndPos,EndOrient,EndType");
            writer.write(",EndGene,EndTranscript,EndRegionType,EndExon,EndPhase,EndCodingBases,EndTotalCodingBases,EndCodingStart,EndCodingEnd");
            writer.newLine();

            for(final GeneFusion fusion : fusions)
            {
                final Transcript startTrans = fusion.upstreamLinkedAnnotation();
                final Transcript endTrans = fusion.downstreamLinkedAnnotation();

                final GeneAnnotation startVar = startTrans.parent();
                final GeneAnnotation endVar = endTrans.parent();

                writer.write(String.format("%s", mSampleId));

                // write upstream SV, transcript and exon info
                writer.write(
                        String.format(",%d,%s,%d,%d,%s",
                                startVar.id(), startVar.chromosome(), startVar.position(), startVar.orientation(), startVar.type()));

                writer.write(
                        String.format(",%s,%s,%s,%d,%d,%d,%d,%d,%d",
                                startTrans.parent().geneName(), startTrans.transcriptId(), startTrans.getRegionType(), startTrans.exonUpstream(), startTrans.exonUpstreamPhase(),
                                startTrans.codingBases(), startTrans.totalCodingBases(),
                                startTrans.codingStart() != null ? startTrans.codingStart() : 0,
                                startTrans.codingEnd() != null ? startTrans.codingEnd() : 0));

                writer.write(
                        String.format(",%d,%s,%d,%d,%s",
                                endVar.id(), endVar.chromosome(), endVar.position(), endVar.orientation(), endVar.type()));

                writer.write(
                        String.format(",%s,%s,%s,%d,%d,%d,%d,%d,%d",
                                endTrans.parent().geneName(), endTrans.transcriptId(), endTrans.getRegionType(), endTrans.exonUpstream(), endTrans.exonUpstreamPhase(),
                                endTrans.codingBases(), endTrans.totalCodingBases(),
                                endTrans.codingStart() != null ? endTrans.codingStart() : 0,
                                endTrans.codingEnd() != null ? endTrans.codingEnd() : 0));

                writer.newLine();
            }

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing gene annotations");
        }
    }


}
