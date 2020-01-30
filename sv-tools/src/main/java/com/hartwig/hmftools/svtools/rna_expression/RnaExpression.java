package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.LinxConfig.formOutputPath;
import static com.hartwig.hmftools.svtools.common.ConfigUtils.DATA_OUTPUT_DIR;
import static com.hartwig.hmftools.svtools.common.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.GENE_TRANSCRIPTS_DIR;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.SAMPLE;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.createCmdLineOptions;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class RnaExpression
{
    private final RnaExpConfig mConfig;
    private final String mSampledId;
    private final RnaBamReader mRnaBamReader;
    private final SvGeneTranscriptCollection mGeneTransCache;
    private final List<GeneReadData> mGeneReadDatalist;
    private BufferedWriter mWriter;

    private static final Logger LOGGER = LogManager.getLogger(RnaExpression.class);

    public RnaExpression(final CommandLine cmd)
    {
        mConfig = new RnaExpConfig(cmd);

        mRnaBamReader = new RnaBamReader(mConfig);

        mSampledId = cmd.getOptionValue(SAMPLE);

        mGeneTransCache = new SvGeneTranscriptCollection();
        mGeneTransCache.setDataPath(cmd.getOptionValue(GENE_TRANSCRIPTS_DIR));


        if(!mConfig.RestrictedGeneIds.isEmpty())
        {
            mGeneTransCache.setRestrictedGeneIdList(mConfig.RestrictedGeneIds);
        }

        mGeneTransCache.setRequiredData(true, false, false, !mConfig.AllTranscripts);
        mGeneTransCache.loadEnsemblData(false);

        mGeneReadDatalist = Lists.newArrayList();
        mWriter = null;
    }

    public void runAnalysis()
    {
        // measure read counts of exonic regions for all specific genes
        int geneCount = 0;
        for(Map.Entry<String,List<EnsemblGeneData>> entry : mGeneTransCache.getChrGeneDataMap().entrySet())
        {
            final List<EnsemblGeneData> genesDataList = entry.getValue();

            for(EnsemblGeneData geneData : genesDataList)
            {
                processGene(geneData);
                ++geneCount;

                if(geneCount > 1 && (geneCount % 100) == 0)
                    LOGGER.info("processed {} genes", geneCount);
            }
        }

        closeBufferedWriter(mWriter);
    }

    private void processGene(final EnsemblGeneData geneData)
    {
        GeneReadData geneReadData = new GeneReadData(geneData);

        final List<TranscriptData> transDataList = mGeneTransCache.getTranscripts(geneData.GeneId);
        if(transDataList.isEmpty())
        {
            LOGGER.warn("no transcripts found for gene({}:{})", geneData.GeneId, geneData.GeneName);
            return;
        }

        // form a genomic region for each unique exon amongst the transcripts
        long minTransPos = -1;
        long maxTransPos = 0;
        RegionReadData prevRegionReadData = null;

        for(final TranscriptData transData : transDataList)
        {
            for(int i = 0; i < transData.exons().size(); ++ i)
            {
                ExonData exon = transData.exons().get(i);

                RegionReadData regionReadData = geneReadData.findRegionData(exon.ExonStart, exon.ExonEnd);

                if (regionReadData == null)
                {
                    GenomeRegion region = GenomeRegions.create(geneData.Chromosome, exon.ExonStart, exon.ExonEnd);
                    regionReadData = new RegionReadData(
                            region, String.format("trans(%s) exon(%d)", transData.TransName, exon.ExonRank));
                    geneReadData.addRegionReadData(regionReadData);
                }

                if(prevRegionReadData != null)
                {
                    prevRegionReadData.getPostRegions().add(regionReadData);
                    regionReadData.getPreRegions().add(prevRegionReadData);
                }

                prevRegionReadData = regionReadData;
            }

            maxTransPos = max(transData.TransEnd, maxTransPos);

            if(minTransPos < 0 || transData.TransStart < minTransPos)
                minTransPos = transData.TransStart;
        }

        GenomeRegion geneRegion = GenomeRegions.create(geneData.Chromosome, minTransPos, maxTransPos);
        mRnaBamReader.readBamCounts(geneReadData, geneRegion);
        mRnaBamReader.analyseReads();

        // report evidence for each gene transcript
        for(final TranscriptData transData : transDataList)
        {
            writeResults(geneReadData, transData);
        }

        // mGeneReadDatalist.add(geneReadData);
    }

    private void writeResults(final GeneReadData geneReadData, final TranscriptData transData)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            if(mWriter == null)
            {
                final String outputFileName = mConfig.OutputDir + "RNA_GENE_EXPRESSION.csv";

                mWriter = createBufferedWriter(outputFileName, false);
                mWriter.write("SampleId,GeneId,GeneName,TransId,ExonCount");
                mWriter.write(",ExonsMatched,AvgDepth,TotalReads,ExonLength,CodingCoverage");
                mWriter.write(",LinksMatched,TotalLinkingReads,NonAdjacentReads");
                mWriter.newLine();
            }

            int exonsFound = 0;
            int linksFound = 0;
            int totalLinkedReads = 0;
            double readDepthTotal = 0;
            long codingLengthTotal = 0;
            long totalReadCoverage = 0;
            int nonAdjacentReads = 0;

            final List<ExonData> exons = transData.exons();

            for(int i = 0; i < exons.size(); ++i)
            {
                ExonData exon = exons.get(i);

                final RegionReadData exonReadData = geneReadData.findRegionData(exon.ExonStart, exon.ExonEnd);
                if(exonReadData == null)
                    continue;

                int exonCoverage = exonReadData.baseCoverage(1);
                totalReadCoverage += exonCoverage;

                if(exonCoverage / (double)exonReadData.length() > 0.50)
                {
                    ++exonsFound;
                }

                codingLengthTotal += exon.ExonEnd - exon.ExonStart;
                nonAdjacentReads += exonReadData.nonAdjacentReads();

                ExonData prevExon = i > 0 ? exons.get(i - 1) : null;
                ExonData nextExon = i < exons.size() - 1 ? exons.get(i + 1) : null;

                RegionReadData prevRegion = prevExon != null ? geneReadData.findRegionData(prevExon.ExonStart, prevExon.ExonEnd) : null;
                RegionReadData nextRegion = nextExon != null ? geneReadData.findRegionData(nextExon.ExonStart, nextExon.ExonEnd) : null;

                boolean linked = (prevExon == null || exonReadData.getLinkedRegions().containsKey(prevRegion))
                        && (nextExon == null || exonReadData.getLinkedRegions().containsKey(nextRegion));

                totalLinkedReads += exonReadData.getLinkedRegions().values().stream().mapToInt(x -> x).sum();

                if(linked)
                    ++linksFound;

                readDepthTotal += exonReadData.averageDepth();
            }

            double avgReadDepth = exonsFound > 0 ? readDepthTotal / exonsFound : 0;

            mWriter.write(String.format("%s,%s,%s,%s,%d",
                    mSampledId, geneReadData.GeneData.GeneId, geneReadData.GeneData.GeneName,
                    transData.TransName, exons.size()));

            mWriter.write(String.format(",%d,%.0f,%d,%d,%d",
                    exonsFound, avgReadDepth, geneReadData.totalReadCount(), codingLengthTotal, totalReadCoverage));

            mWriter.write(String.format(",%d,%d,%d",
                    linksFound, totalLinkedReads, nonAdjacentReads));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            LOGGER.error("failed to write gene expression file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = createCmdLineOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        if(!RnaExpConfig.checkValid(cmd))
        {
            LOGGER.error("missing config options, exiting");
            return;
        }

        RnaExpression rnaExpression = new RnaExpression(cmd);
        rnaExpression.runAnalysis();

        LOGGER.info("RNA expression analysis complete");
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
