package com.hartwig.hmftools.linx.visualiser.file;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.analysis.SvClassification.isFilteredResolvedType;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_Q;
import static com.hartwig.hmftools.linx.types.SvBreakend.DIRECTION_CENTROMERE;
import static com.hartwig.hmftools.linx.types.SvBreakend.DIRECTION_TELOMERE;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.isStart;
import static com.hartwig.hmftools.linx.visualiser.file.VisProteinDomainFile.PD_FIVE_PRIME_UTR;
import static com.hartwig.hmftools.linx.visualiser.file.VisProteinDomainFile.PD_NON_CODING;
import static com.hartwig.hmftools.linx.visualiser.file.VisProteinDomainFile.PD_THREE_PRIME_UTR;
import static com.hartwig.hmftools.linx.visualiser.file.VisSvDataFile.INFO_TYPE_FOLDBACK;
import static com.hartwig.hmftools.linx.visualiser.file.VisSvDataFile.INFO_TYPE_NORMAL;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptProteinData;
import com.hartwig.hmftools.linx.cn.SvCNData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.chaining.SvChain;
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

    private boolean mEnabled;
    private final String mOutputDir;
    private String mSampleId;

    private boolean mBatchOutput;
    private BufferedWriter mSvFileWriter;
    private BufferedWriter mSegmentFileWriter;
    private BufferedWriter mCnFileWriter;
    private BufferedWriter mGeneFileWriter;
    private BufferedWriter mProteinDomainFileWriter;

    private static final Logger LOGGER = LogManager.getLogger(VisualiserWriter.class);


    public VisualiserWriter(final String outputDir, boolean enabled, boolean isBatchOutput)
    {
        mEnabled = enabled;
        mOutputDir = outputDir;
        mSampleId = "";
        mBatchOutput = isBatchOutput;

        if(mBatchOutput && mEnabled)
        {
            initialiseBatchOutputFiles();
        }

        mGeneTranscriptCollection = null;

        mGeneData = Lists.newArrayList();
    }

    private void initialiseBatchOutputFiles()
    {
        try
        {
            mSvFileWriter = createBufferedWriter(mOutputDir + "SVA_VIS_SVS.tsv", false);
            mSvFileWriter.write(VisSvDataFile.header());
            mSvFileWriter.newLine();

            mSegmentFileWriter = createBufferedWriter(mOutputDir + "SVA_VIS_SEGMENTS.tsv", false);
            mSegmentFileWriter.write(VisSegmentFile.header());
            mSegmentFileWriter.newLine();

            mCnFileWriter = createBufferedWriter(mOutputDir + "SVA_VIS_COPY_NUMBER.tsv", false);
            mCnFileWriter.write(VisCopyNumberFile.header());
            mCnFileWriter.newLine();

            mGeneFileWriter = createBufferedWriter(mOutputDir + "SVA_VIS_GENE_EXONS.tsv", false);
            mGeneFileWriter.write(VisGeneExonFile.header());
            mGeneFileWriter.newLine();

            mProteinDomainFileWriter = createBufferedWriter(mOutputDir + "SVA_VIS_PROTEIN_DOMAINS.tsv", false);
            mProteinDomainFileWriter.write(VisProteinDomainFile.header());
            mProteinDomainFileWriter.newLine();
        }
        catch(IOException e)
        {
            LOGGER.error("failed to open and write output file headers");
        }
    }

    public void close()
    {
        closeBufferedWriter(mSvFileWriter);
        closeBufferedWriter(mSegmentFileWriter);
        closeBufferedWriter(mCnFileWriter);
        closeBufferedWriter(mGeneFileWriter);
        closeBufferedWriter(mProteinDomainFileWriter);
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
        writeGeneData();
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

            // int repeatCount = (int)var.getRoundedPloidy(true);

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
                        var.ploidy()));
            }
        }

        try
        {
            if(mBatchOutput)
            {
                for(final VisSvDataFile data : svDataList)
                {
                    mSvFileWriter.write(VisSvDataFile.toString(data));
                    mSvFileWriter.newLine();
                }
            }
            else
            {
                VisSvDataFile.write(VisSvDataFile.generateFilename(mOutputDir, mSampleId), svDataList);
            }
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
                boolean startsOnEnd = chain.getFirstSV().equals(chain.getLastSV(), true); // closed loop chains eg DMs
                double chainPloidy = chain.ploidy();

                if(!chain.isClosedLoop())
                {
                    SvBreakend breakend = chain.getOpenBreakend(true);

                    if (breakend != null)
                    {
                        segments.add(new VisSegmentFile(mSampleId, cluster.id(), chain.id(), breakend.chromosome(),
                                getPositionValue(breakend, true), getPositionValue(breakend, false), chainPloidy));
                    }
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

                    if(chainPloidy == 0)
                    {
                        int pairRepeatCount = pair.repeatCount();

                        // check for duplicate links - would only exist in non-identical chains since these have already been removed
                        for (final SvChain otherChain : cluster.getChains())
                        {
                            pairRepeatCount += otherChain.getLinkedPairs().stream().filter(x -> x.matches(pair)).count();
                        }

                        chainPloidy = pairRepeatCount;
                    }

                    final SvBreakend beStart = pair.getBreakend(true);
                    final SvBreakend beEnd = pair.getBreakend(false);

                    segments.add(new VisSegmentFile(mSampleId, cluster.id(), chain.id(),
                            beStart.chromosome(), Long.toString(beStart.position()), Long.toString(beEnd.position()), chainPloidy));
                }

                if(!chain.isClosedLoop())
                {
                    // log the end of the chain out to centromere or telomere
                    SvBreakend breakend = chain.getOpenBreakend(false);

                    if (breakend != null && !startsOnEnd)
                    {
                        segments.add(new VisSegmentFile(mSampleId, cluster.id(), chain.id(), breakend.chromosome(),
                                getPositionValue(breakend, true), getPositionValue(breakend, false), chainPloidy));
                    }
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
            if(mBatchOutput)
            {
                for(final VisSegmentFile data : segments)
                {
                    mSegmentFileWriter.write(VisSegmentFile.toString(data));
                    mSegmentFileWriter.newLine();
                }
            }
            else
            {
                VisSegmentFile.write(VisSegmentFile.generateFilename(mOutputDir, mSampleId), segments);
            }
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


    public void addGeneExonData(int clusterId, final String geneId, final String geneName, final String transName, int transId,
            final String chromosome, final String annotationType)
    {
        mGeneData.add(new VisGeneData(clusterId, geneId, geneName, transName, transId, chromosome, annotationType));
    }

    private void writeGeneData()
    {
        if(!mEnabled)
            return;

        List<VisGeneExonFile> geneExonList = Lists.newArrayList();
        List<VisProteinDomainFile> proteinList = Lists.newArrayList();

        // exons: SampleId,ClusterId,Gene,Transcript,Chromosome,AnnotationType,ExonRank,ExonStart,ExonEnd

        // protein domains: SampleId,ClusterId,Gene,Transcript,Start,End,Info

        // first remove duplicates from amongst the genes
        List<String> loggedTranscripts = Lists.newArrayList();

        for(final VisGeneData geneData : mGeneData)
        {
            if(geneData.ClusterId < 0) // skip drivers not linked to a cluster
                continue;

            final TranscriptData transData = mGeneTranscriptCollection.getTranscriptData(geneData.GeneId, geneData.TransName);

            if (loggedTranscripts.contains(transData.TransName))
                continue;

            loggedTranscripts.add(transData.TransName);

            if(transData == null || transData.exons().isEmpty())
                continue;

            for (final ExonData exonData : transData.exons())
            {
                geneExonList.add(new VisGeneExonFile(mSampleId, geneData.ClusterId, geneData.GeneName, transData.TransName,
                        geneData.Chromosome, geneData.AnnotationType, exonData.ExonRank, exonData.ExonStart, exonData.ExonEnd));
            }

            int transId = geneData.TransId > 0 ? geneData.TransId : transData.TransId;

            final List<TranscriptProteinData> transProteinData = mGeneTranscriptCollection.getTranscriptProteinDataMap().get(transId);

            if (transProteinData != null)
            {
                for (final TranscriptProteinData proteinData : transProteinData)
                {
                    final Long[] domainPositions = mGeneTranscriptCollection.getProteinDomainPositions(proteinData, transData);

                    if(domainPositions[SE_START] != null && domainPositions[SE_END] != null)
                    {
                        proteinList.add(new VisProteinDomainFile(mSampleId, geneData.ClusterId, transData.TransName, geneData.Chromosome,
                                domainPositions[SE_START], domainPositions[SE_END], proteinData.HitDescription));
                    }
                }
            }

            // show the 5' and 3' UTR or non-coding regions as 'protein domains'
            if(transData.CodingEnd != null && transData.CodingEnd != null)
            {
                long fivePrimeUtrStart = transData.Strand == 1 ? transData.TransStart : transData.CodingEnd + 1;
                long fivePrimeUtrEnd = transData.Strand == 1 ? transData.CodingStart - 1 : transData.TransEnd;

                long threePrimeUtrStart = transData.Strand == 1 ? transData.CodingEnd + 1 : transData.TransStart;
                long threePrimeUtrEnd = transData.Strand == 1 ? transData.TransEnd : transData.CodingStart - 1;

                if(fivePrimeUtrStart < fivePrimeUtrEnd)
                {
                    proteinList.add(new VisProteinDomainFile(mSampleId, geneData.ClusterId, transData.TransName, geneData.Chromosome,
                            fivePrimeUtrStart, fivePrimeUtrEnd, PD_FIVE_PRIME_UTR));
                }

                if(threePrimeUtrStart < threePrimeUtrEnd)
                {
                    proteinList.add(new VisProteinDomainFile(mSampleId, geneData.ClusterId, transData.TransName, geneData.Chromosome,
                            threePrimeUtrStart, threePrimeUtrEnd, PD_THREE_PRIME_UTR));
                }
            }
            else
            {
                proteinList.add(new VisProteinDomainFile(mSampleId, geneData.ClusterId, transData.TransName, geneData.Chromosome,
                        transData.TransStart, transData.TransEnd, PD_NON_CODING));
            }
        }

        try
        {
            if(mBatchOutput)
            {
                for(final VisGeneExonFile data : geneExonList)
                {
                    mGeneFileWriter.write(VisGeneExonFile.toString(data));
                    mGeneFileWriter.newLine();
                }

                for(final VisProteinDomainFile data : proteinList)
                {
                    mProteinDomainFileWriter.write(VisProteinDomainFile.toString(data));
                    mProteinDomainFileWriter.newLine();
                }
            }
            else
            {
                VisGeneExonFile.write(VisGeneExonFile.generateFilename(mOutputDir, mSampleId), geneExonList);
                VisProteinDomainFile.write(VisProteinDomainFile.generateFilename(mOutputDir, mSampleId), proteinList);
            }
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

        try
        {
            if(mBatchOutput)
            {
                for(final VisCopyNumberFile data : cnDataList)
                {
                    mCnFileWriter.write(VisCopyNumberFile.toString(data));
                    mCnFileWriter.newLine();
                }
            }
            else
            {
                VisCopyNumberFile.write(VisCopyNumberFile.generateFilename(mOutputDir, mSampleId), cnDataList);
            }
        }
        catch(IOException e)
        {
            LOGGER.error("filed to write VIS copy number output");
        }
    }
}
