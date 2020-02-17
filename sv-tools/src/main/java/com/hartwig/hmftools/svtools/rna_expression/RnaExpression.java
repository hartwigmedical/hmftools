package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.common.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TC_UNSPLICED;
import static com.hartwig.hmftools.svtools.rna_expression.RegionReadData.findUniqueBases;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.GENE_FRAGMENT_BUFFER;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.GENE_TRANSCRIPTS_DIR;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.SAMPLE;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.createCmdLineOptions;
import static com.hartwig.hmftools.svtools.rna_expression.TranscriptModel.calculateTranscriptResults;
import static com.hartwig.hmftools.svtools.rna_expression.TranscriptModel.estimateRatesByLeastSquares;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
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

public class RnaExpression
{
    private final RnaExpConfig mConfig;
    private final String mSampledId;
    private final RnaBamReader mRnaBamReader;
    private final ResultsWriter mResultsWriter;
    private final SvGeneTranscriptCollection mGeneTransCache;
    private final GcBiasAdjuster mGcBiasAdjuster;
    private final FragmentSizeCalcs mFragmentSizeCalcs;
    private final ExpectedExpressionRates mExpExpressionRates;

    private static final Logger LOGGER = LogManager.getLogger(RnaExpression.class);

    public RnaExpression(final CommandLine cmd)
    {
        mConfig = new RnaExpConfig(cmd);

        mRnaBamReader = new RnaBamReader(mConfig);
        mGcBiasAdjuster = new GcBiasAdjuster(mConfig);

        mSampledId = cmd.getOptionValue(SAMPLE);
        mResultsWriter = new ResultsWriter(mConfig);
        mResultsWriter.setSampleId(mSampledId);

        mGeneTransCache = new SvGeneTranscriptCollection();
        mGeneTransCache.setDataPath(cmd.getOptionValue(GENE_TRANSCRIPTS_DIR));

        if(!mConfig.RestrictedGeneIds.isEmpty())
        {
            mGeneTransCache.setRestrictedGeneIdList(mConfig.RestrictedGeneIds);
        }

        mGeneTransCache.setRequiredData(true, false, false, mConfig.CanonicalTranscriptOnly);
        mGeneTransCache.loadEnsemblData(false);

        mFragmentSizeCalcs = new FragmentSizeCalcs(mConfig, mGeneTransCache, mRnaBamReader);

        mExpExpressionRates = mConfig.GenerateExpectedExpression ? new ExpectedExpressionRates(mConfig) : null;
    }

    public void runAnalysis()
    {
        LOGGER.info("sample({}) running RNA expression analysis", mSampledId);

        if(!mRnaBamReader.validReader())
        {
            LOGGER.warn("BAM reader init failed");
            return;
        }

        if(mGcBiasAdjuster.enabled())
        {
            mGcBiasAdjuster.loadData();
            mGcBiasAdjuster.generateDepthCounts(mRnaBamReader, mGeneTransCache.getChrGeneDataMap());
        }

        if(mConfig.FragmentLengthMinCount > 0)
        {
            mFragmentSizeCalcs.calcSampleFragmentSize();

            if (mConfig.WriteFragmentLengths)
                return;
        }

        // measure read counts of exonic regions for all specific genes
        int geneCount = 0;
        for(Map.Entry<String,List<EnsemblGeneData>> entry : mGeneTransCache.getChrGeneDataMap().entrySet())
        {
            final List<EnsemblGeneData> genesDataList = entry.getValue();

            for(EnsemblGeneData geneData : genesDataList)
            {
                if(mConfig.ExcludedGeneIds.contains(geneData.GeneId))
                    continue;

                processGene(geneData);
                ++geneCount;

                if(geneCount > 1 && (geneCount % 100) == 0)
                    LOGGER.info("processed {} genes", geneCount);
            }
        }

        mResultsWriter.close();
        mRnaBamReader.close();

        if(mExpExpressionRates != null)
            mExpExpressionRates.close();
    }

    private void processGene(final EnsemblGeneData geneData)
    {
        GeneReadData geneReadData = new GeneReadData(geneData);

        List<TranscriptData> transDataList = Lists.newArrayList(mGeneTransCache.getTranscripts(geneData.GeneId));

        if(transDataList.isEmpty())
        {
            LOGGER.warn("no transcripts found for gene({}:{})", geneData.GeneId, geneData.GeneName);
            return;
        }

        if(!mConfig.SpecificTransIds.isEmpty())
            transDataList = transDataList.stream().filter(x -> mConfig.SpecificTransIds.contains(x.TransName)).collect(Collectors.toList());

        geneReadData.setTranscripts(transDataList);

        geneReadData.generateExonicRegions();

        // cache reference bases for comparison with read bases
        if(mConfig.RefFastaSeqFile != null)
        {
            for (RegionReadData region : geneReadData.getExonRegions())
            {
                final String regionRefBases = mConfig.RefFastaSeqFile.getSubsequenceAt(
                        region.chromosome(), region.start(), region.end()).getBaseString();

                region.setRefBases(regionRefBases);
            }

            findUniqueBases(geneReadData.getExonRegions());
        }

        // use a buffer around the gene to pick up reads which span outside its transcripts
        long regionStart = geneReadData.getTranscriptsRange()[SE_START] - GENE_FRAGMENT_BUFFER;
        long regionEnd = geneReadData.getTranscriptsRange()[SE_END] + GENE_FRAGMENT_BUFFER;

        if(regionStart >= regionEnd)
        {
            LOGGER.warn("invalid gene({}:{}) region({} -> {})", geneData.GeneId, geneData.GeneName, regionStart, regionEnd);
            return;
        }

        GenomeRegion geneRegion = GenomeRegions.create(geneData.Chromosome, regionStart, regionEnd);

        mRnaBamReader.readBamCounts(geneReadData, geneRegion);

        runTranscriptEstimation(geneReadData);

        mResultsWriter.writeGeneData(geneReadData);

        if(!mConfig.GeneStatsOnly)
        {
            // report evidence for each gene transcript
            for (final TranscriptData transData : transDataList)
            {
                final TranscriptResults results = calculateTranscriptResults(geneReadData, transData);
                geneReadData.getTranscriptResults().add(results);

                mResultsWriter.writeTranscriptResults(geneReadData, results);

                if (mConfig.WriteExonData)
                {
                    mResultsWriter.writeExonData(geneReadData, transData);
                }
            }
        }
    }

    private void runTranscriptEstimation(final GeneReadData geneReadData)
    {
        if(mExpExpressionRates == null)
            return;

        mExpExpressionRates.generateExpectedRates(geneReadData);

        final double[] transcriptCounts = mExpExpressionRates.generateTranscriptCounts(
                geneReadData, mRnaBamReader.getTransComboData(), geneReadData.getCounts()[TC_UNSPLICED]);

        estimateRatesByLeastSquares(transcriptCounts, mExpExpressionRates.getTranscriptDefinitions(), mExpExpressionRates.getTranscriptNames());
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
