package com.hartwig.hmftools.linx.visualiser.file;

import static java.lang.Math.max;

import static com.hartwig.hmftools.linx.analysis.SvClassification.isFilteredResolvedType;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_Q;
import static com.hartwig.hmftools.linx.types.SvBreakend.DIRECTION_CENTROMERE;
import static com.hartwig.hmftools.linx.types.SvBreakend.DIRECTION_TELOMERE;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.isStart;
import static com.hartwig.hmftools.linx.visualiser.file.VisSvDataFile.INFO_TYPE_FOLDBACK;
import static com.hartwig.hmftools.linx.visualiser.file.VisSvDataFile.INFO_TYPE_NORMAL;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptExonData;
import com.hartwig.hmftools.linx.cn.SvCNData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvChain;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class VisualiserWriter
{
    // references only
    SvGeneTranscriptCollection mGeneTranscriptCollection;

    private List<VisGeneData> mGeneData;

    private boolean mBatchOutput;
    private boolean mEnabled;
    private final String mOutputDir;
    private String mSampleId;

    private static final Logger LOGGER = LogManager.getLogger(VisualiserWriter.class);


    public VisualiserWriter(final String outputDir, boolean enabled, boolean isBatchOutput)
    {
        mEnabled = enabled;
        mOutputDir = outputDir;
        mSampleId = "";
        mBatchOutput = isBatchOutput;

        mGeneTranscriptCollection = null;

        mGeneData = Lists.newArrayList();
    }

    public void setGeneDataCollection(SvGeneTranscriptCollection geneTranscriptCollection)
    {
        mGeneTranscriptCollection = geneTranscriptCollection;
    }

    public void setSampleId(final String sampleId)
    {
        mSampleId = sampleId;
        mGeneData.clear();
    }

    public void writeOutput(final List<SvCluster> clusters, final List<SvVarData> variants, final Map<String,List<SvCNData>> chrCnDataMap)
    {
        if(!mEnabled)
            return;

        writeVisualSvData(variants);
        writeVisualSegmentData(clusters);
        writeGeneExonData();
        writeCopyNumberData(chrCnDataMap);
    }

    private void writeVisualSvData(final List<SvVarData> variants)
    {
        List<VisSvDataFile> svDataList = Lists.newArrayList();

        for(final SvVarData var : variants)
        {
            final List<SvChain> chains = var.getCluster().findChains(var);

            // repeat an SV for every time it appears in a chain
            int chainCount = chains.isEmpty() ? 1 : chains.size();
            int unchainedChainId = chains.isEmpty() ? var.getCluster().getChainId(var) : -1;

            int repeatCount = max((int)var.getRoundedCNChange(),1);

            for(int i = 0; i < chainCount; ++i)
            {
                int chainId = chains.isEmpty() ? unchainedChainId : chains.get(i).id();

                final SvBreakend beStart = var.getBreakend(true);
                final SvBreakend beEnd = var.getBreakend(false);

                svDataList.add(new VisSvDataFile(mSampleId, var.getCluster().id(), chainId, var.dbId(),
                        var.type(), var.getCluster().getResolvedType(),
                        beStart.chromosome(), beEnd != null ? beEnd.chromosome() : "-1",
                        beStart.position(),beEnd != null ? beEnd.position() : 0,
                        beStart.orientation(), beEnd != null ? beEnd.orientation() : 0,
                        beStart.getSV().getFoldbackLink(beStart.usesStart()).isEmpty() ? INFO_TYPE_NORMAL : INFO_TYPE_FOLDBACK,
                        beEnd!= null ? (beEnd.getSV().getFoldbackLink(beEnd.usesStart()).isEmpty() ? INFO_TYPE_NORMAL : INFO_TYPE_FOLDBACK) : "",
                        repeatCount));
            }
        }

        final String outputFile = mBatchOutput ? mOutputDir + "SVA_VIS_SVS.csv" : VisSvDataFile.generateFilename(mOutputDir, mSampleId);

        try
        {
            VisSvDataFile.write(outputFile, svDataList);
        }
        catch(IOException e)
        {
            LOGGER.error("filed to write VIS copy number output");
        }
    }

    private void writeVisualSegmentData(final List<SvCluster> clusters)
    {
        List<VisSegmentFile> segments = Lists.newArrayList();

        for(final SvCluster cluster : clusters)
        {
            if (cluster.getSvCount() == 1 || isFilteredResolvedType(cluster.getResolvedType()))
                continue;

            // for any linked pair which is repeated in a separate chain, skip writing it for subsequent chains
            List<SvLinkedPair> uniquePairs = Lists.newArrayList();

            for (final SvChain chain : cluster.getChains())
            {
                // log the start of the chain
                SvBreakend breakend = chain.getOpenBreakend(true);
                boolean startsOnEnd = chain.getFirstSV().equals(chain.getLastSV(), true);

                if (breakend != null)
                {
                    segments.add(new VisSegmentFile(mSampleId, cluster.id(), chain.id(), breakend.chromosome(),
                            getPositionValue(breakend, true), getPositionValue(breakend, false),
                            startsOnEnd ? 2 : 1));
                }

                for (final SvLinkedPair pair : chain.getLinkedPairs())
                {
                    boolean isRepeat = false;

                    // only log each chain link once, and log how many times the link has been used
                    for (final SvLinkedPair existingPair : uniquePairs)
                    {
                        if (pair.matches(existingPair))
                        {
                            isRepeat = true;
                            break;
                        }
                    }

                    if (isRepeat)
                        continue;

                    uniquePairs.add(pair);

                    int pairRepeatCount = 0;

                    for (final SvChain otherChain : cluster.getChains())
                    {
                        pairRepeatCount += otherChain.getLinkedPairs().stream().filter(x -> x.matches(pair)).count();
                    }

                    final SvBreakend beStart = pair.getBreakend(true);
                    final SvBreakend beEnd = pair.getBreakend(false);

                    segments.add(new VisSegmentFile(mSampleId, cluster.id(), chain.id(),
                            beStart.chromosome(), Long.toString(beStart.position()), Long.toString(beEnd.position()), pairRepeatCount));
                }

                // log the end of the chain out to centromere or telomere
                breakend = chain.getOpenBreakend(false);

                if (breakend != null && !startsOnEnd)
                {
                    segments.add(new VisSegmentFile(mSampleId, cluster.id(), chain.id(), breakend.chromosome(),
                            getPositionValue(breakend, true), getPositionValue(breakend, false), 1));
                }
            }

            // finally write out all unchained SVs
            for (final SvVarData var : cluster.getUnlinkedSVs())
            {
                if (var.isReplicatedSv())
                    continue;

                int chainId = cluster.getChainId(var);

                for (int be = SE_START; be <= SE_END; ++be)
                {
                    final SvBreakend breakend = var.getBreakend(isStart(be));

                    if (breakend == null)
                        continue;

                    segments.add(new VisSegmentFile(mSampleId, cluster.id(), chainId, breakend.chromosome(),
                            getPositionValue(breakend, true), getPositionValue(breakend, false), 1));
                }
            }
        }

        try
        {
            String outputFileName = mBatchOutput ? mOutputDir + "SVA_VIS_SEGMENTS.csv" :
                    VisSegmentFile.generateFilename(mOutputDir, mSampleId);

            VisSegmentFile.write(outputFileName, segments);
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to visualiser segments file: {}", e.toString());
        }
    }

    private static final String getPositionValue(final SvBreakend breakend, boolean isChainEnd)
    {
        if(breakend.orientation() == 1 && breakend.arm().equals(CHROMOSOME_ARM_P))
        {
            return isChainEnd ? DIRECTION_TELOMERE : Long.toString(breakend.position());
        }
        else if(breakend.orientation() == -1 && breakend.arm().equals(CHROMOSOME_ARM_P))
        {
            return isChainEnd ? Long.toString(breakend.position()) : DIRECTION_CENTROMERE;
        }
        else if(breakend.orientation() == -1 && breakend.arm().equals(CHROMOSOME_ARM_Q))
        {
            return isChainEnd ? Long.toString(breakend.position()) : DIRECTION_TELOMERE;
        }
        else
        {
            return isChainEnd ? DIRECTION_CENTROMERE : Long.toString(breakend.position());
        }
    }


    public void addGeneExonData(int clusterId, final String geneId, final String geneName, final String transcriptId,
            final String chromosome, final String annotationType)
    {
        mGeneData.add(new VisGeneData(clusterId, geneId, geneName, transcriptId, chromosome, annotationType));
    }

    public void writeGeneExonData()
    {
        if(!mEnabled)
            return;

        List<VisGeneExonFile> geneExonList = Lists.newArrayList();

        // SampleId,ClusterId,Gene,Transcript,Chromosome,AnnotationType,ExonRank,ExonStart,ExonEnd

        // first remove duplicates from amongst the genes
        List<String> loggedGenes = Lists.newArrayList();

        for(final VisGeneData geneData : mGeneData)
        {
            if (loggedGenes.contains(geneData.GeneId))
                continue;

            loggedGenes.add(geneData.GeneId);

            final List<TranscriptExonData> exonDataLst =
                    mGeneTranscriptCollection.getTranscriptExons(geneData.GeneId, geneData.TranscriptId);

            for (final TranscriptExonData exonData : exonDataLst)
            {
                geneExonList.add(new VisGeneExonFile(mSampleId, geneData.ClusterId, geneData.GeneName, exonData.TransName,
                        geneData.Chromosome, geneData.AnnotationType, exonData.ExonRank, exonData.ExonStart, exonData.ExonEnd));
            }
        }

        try
        {
            String outputFileName = mBatchOutput ? mOutputDir + "SVA_VIS_GENE_EXONS.csv" :
                    VisGeneExonFile.generateFilename(mOutputDir, mSampleId);

            VisGeneExonFile.write(outputFileName, geneExonList);

        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to visualiser gene-exons file: {}", e.toString());
        }
    }

    private void writeCopyNumberData(final Map<String,List<SvCNData>> chrCNDataMap)
    {
        List<VisCopyNumberFile> cnDataList = Lists.newArrayList();

        for(Map.Entry<String,List<SvCNData>> entry : chrCNDataMap.entrySet())
        {
            final String chromosome = entry.getKey();

            for(SvCNData cnData : entry.getValue())
            {
                cnDataList.add(new VisCopyNumberFile(
                        mSampleId, chromosome, cnData.StartPos, cnData.EndPos, cnData.CopyNumber, cnData.ActualBaf));
            }
        }

        final String outputFile = mBatchOutput ? mOutputDir + "SVA_VIS_COPY_NUMBER.csv" :
                VisCopyNumberFile.generateFilename(mOutputDir, mSampleId);

        try
        {
            VisCopyNumberFile.write(outputFile, cnDataList);
        }
        catch(IOException e)
        {
            LOGGER.error("filed to write VIS copy number output");
        }
    }
}
