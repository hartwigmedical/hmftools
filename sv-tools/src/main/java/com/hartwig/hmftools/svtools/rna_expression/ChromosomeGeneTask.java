package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStr;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.RESIDUAL_PERC;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.RESIDUAL_TOTAL;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.calcResiduals;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.calculateFittedCounts;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.sumVector;
import static com.hartwig.hmftools.svtools.rna_expression.ExpectedExpressionRates.UNSPLICED_CAT_INDEX;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.typeAsInt;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.markOverlappingGeneRegions;
import static com.hartwig.hmftools.svtools.rna_expression.RegionReadData.findUniqueBases;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.GENE_FRAGMENT_BUFFER;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.RE_LOGGER;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpUtils.positionsOverlap;
import static com.hartwig.hmftools.svtools.rna_expression.TranscriptModel.calculateTranscriptResults;

import java.io.BufferedWriter;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletableFuture;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;

public class ChromosomeGeneTask implements Callable
{
    private final String mChromosome;
    private final RnaExpConfig mConfig;
    private final SvGeneTranscriptCollection mGeneTransCache;
    private final ResultsWriter mResultsWriter;

    private final RnaBamReader mRnaBamReader;
    private final ExpectedExpressionRates mExpExpressionRates;
    private final List<EnsemblGeneData> mGeneDataList;
    private int mCurrentGeneIndex;
    public int mGenesProcessed;

    private final PerformanceCounter mPerfCounter;

    public ChromosomeGeneTask(
            final RnaExpConfig config, final String chromosome, final List<EnsemblGeneData> geneDataList,
            final SvGeneTranscriptCollection geneTransCache, final ResultsWriter resultsWriter)
    {
        mConfig = config;
        mChromosome = chromosome;
        mGeneTransCache = geneTransCache;
        mResultsWriter = resultsWriter;

        mGeneDataList = geneDataList;

        mCurrentGeneIndex = 0;
        mRnaBamReader = new RnaBamReader(mConfig);
        mExpExpressionRates = mConfig.GenerateExpectedExpression ? new ExpectedExpressionRates(mConfig) : null;

        mPerfCounter = new PerformanceCounter(String.format("chr(%s) genes(%d)", mChromosome, mGeneDataList.size()));
    }

    public void setWriters(BufferedWriter expRatesWriter, BufferedWriter readsWriter)
    {
        if(mExpExpressionRates == null)
            return;

        mExpExpressionRates.setWriter(expRatesWriter);
        // mRnaBamReader.setWriter(readsWriter);
    }

    public final RnaBamReader getBamReader() { return mRnaBamReader; }

    @Override
    public Long call()
    {
        analyseGenes();
        return (long)1; // return value not used
    }

    public void analyseGenes()
    {
        if(mConfig.RestrictedGeneIds.isEmpty())
        {
            RE_LOGGER.info("processing {} genes for chromosome({})", mGeneDataList.size(), mChromosome);
        }

        while(mCurrentGeneIndex < mGeneDataList.size())
        {
            final List<EnsemblGeneData> overlappingGenes = findNextOverlappingGenes();
            final List<GeneReadData> geneReadDataList = createGeneReadData(overlappingGenes);

            for (GeneReadData geneReadData : geneReadDataList)
            {
                mPerfCounter.start();
                processGene(geneReadData);
                mPerfCounter.stop();

                ++mGenesProcessed;

                if (mGenesProcessed > 1 && (mGenesProcessed % 100) == 0)
                {
                    RE_LOGGER.info("chr({}) processed {} genes", mChromosome, mGenesProcessed);
                }
            }
        }
    }

    public void logPerfStats()
    {
        mPerfCounter.logStats();
    }

    private List<EnsemblGeneData> findNextOverlappingGenes()
    {
        final List<EnsemblGeneData> overlappingGenes = Lists.newArrayList();

        while(mCurrentGeneIndex < mGeneDataList.size())
        {
            EnsemblGeneData geneData = mGeneDataList.get(mCurrentGeneIndex);

            if(mConfig.ExcludedGeneIds.contains(geneData.GeneId))
            {
                ++mCurrentGeneIndex;
                continue;
            }

            if(overlappingGenes.isEmpty() || overlappingGenes.stream().anyMatch(x -> positionsOverlap(geneData.GeneStart, geneData.GeneEnd, x.GeneStart, x.GeneEnd)))
            {
                overlappingGenes.add(geneData);
                ++mCurrentGeneIndex;
            }
            else
            {
                break;
            }
        }

        return overlappingGenes;
    }

    private List<GeneReadData> createGeneReadData(final List<EnsemblGeneData> geneDataList)
    {
        List<GeneReadData> geneReadDataList = Lists.newArrayList();

        for(EnsemblGeneData geneData : geneDataList)
        {
            GeneReadData geneReadData = new GeneReadData(geneData);

            List<TranscriptData> transDataList = Lists.newArrayList(mGeneTransCache.getTranscripts(geneData.GeneId));

            if(transDataList.isEmpty())
            {
                RE_LOGGER.warn("no transcripts found for gene({}:{})", geneData.GeneId, geneData.GeneName);
                continue;
            }

            if(!mConfig.SpecificTransIds.isEmpty())
                transDataList = transDataList.stream().filter(x -> mConfig.SpecificTransIds.contains(x.TransName)).collect(Collectors.toList());

            geneReadData.setTranscripts(transDataList);

            geneReadData.generateExonicRegions();

            geneReadDataList.add(geneReadData);
        }

        if(mConfig.RestrictedGeneIds.isEmpty())
        {
            markOverlappingGeneRegions(geneReadDataList, false);
        }

        // if(geneReadDataList.size() > 1)
        //    logOverlappingGenes(geneReadDataList);

        return geneReadDataList;
    }

    private void logOverlappingGenes(final List<GeneReadData> overlappingGenes)
    {
        String geneNamesStr = "";
        int transcriptCount = 0;
        long minRange = -1;
        long maxRange = 0;

        for(GeneReadData geneReadData : overlappingGenes)
        {
            geneNamesStr = appendStr(geneNamesStr, geneReadData.GeneData.GeneName, ';');
            transcriptCount += geneReadData.getTranscripts().size();
            maxRange =  max(maxRange, geneReadData.GeneData.GeneEnd);
            minRange =  minRange < 0 ? geneReadData.GeneData.GeneStart : min(minRange, geneReadData.GeneData.GeneStart);
        }

        // Time,Chromosome,GeneCount,TranscriptCount,RangeStart,RangeEnd,GeneNames
        RE_LOGGER.info("GENE_OVERLAP: {},{},{},{},{},{}", // chr({}) genes({}) transcripts({}) range({} -> {}),
                mChromosome, overlappingGenes.size(), transcriptCount, minRange, maxRange, geneNamesStr);
    }

    private void processGene(final GeneReadData geneReadData)
    {
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

        final EnsemblGeneData geneData = geneReadData.GeneData;

        if(regionStart >= regionEnd)
        {
            RE_LOGGER.warn("invalid gene({}:{}) region({} -> {})", geneData.GeneId, geneData.GeneName, regionStart, regionEnd);
            return;
        }

        GenomeRegion geneRegion = GenomeRegions.create(geneData.Chromosome, regionStart, regionEnd);

        mRnaBamReader.readBamCounts(geneReadData, geneRegion);

        runTranscriptEstimation(geneReadData);

        mResultsWriter.writeGeneData(geneReadData);

        if(!mConfig.GeneStatsOnly)
        {
            // report evidence for each gene transcript
            for (final TranscriptData transData : geneReadData.getTranscripts())
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

        if(!mExpExpressionRates.validData())
        {
            RE_LOGGER.debug("gene({}) invalid expected rates or actuals data", geneReadData.name());
            return;
        }

        final double[] transComboCounts = mExpExpressionRates.generateTranscriptCounts(geneReadData, mRnaBamReader.getTransComboData());

        double totalCounts = sumVector(transComboCounts);

        if(totalCounts == 0)
            return;

        // add in counts for the unspliced category
        int unsplicedCount = geneReadData.getCounts()[typeAsInt(GeneMatchType.UNSPLICED)];
        transComboCounts[UNSPLICED_CAT_INDEX] = unsplicedCount;

        final List<String> transcriptNames = mExpExpressionRates.getTranscriptNames();

        /*
        final double[] lsqFitAllocations = allocateTranscriptCountsByLeastSquares(transComboCounts, mExpExpressionRates.getTranscriptDefinitions());
        final double[] lsqFittedCounts = calculateFittedCounts(mExpExpressionRates.getTranscriptDefinitions(), lsqFitAllocations);
        double[] lsqResiduals = calcResiduals(transComboCounts, lsqFittedCounts, totalCounts);
        */

        final double[] fitAllocations = ExpectationMaxFit.performFit(transComboCounts, mExpExpressionRates.getTranscriptDefinitions());
        final double[] fittedCounts = calculateFittedCounts(mExpExpressionRates.getTranscriptDefinitions(), fitAllocations);

        double[] residuals = calcResiduals(transComboCounts, fittedCounts, totalCounts);

        RE_LOGGER.debug(String.format("gene(%s) totalFragments(%.0f) residuals(%.0f perc=%.3f)",
                geneReadData.name(), totalCounts, residuals[RESIDUAL_TOTAL], residuals[RESIDUAL_PERC]));

        geneReadData.setFitResiduals(residuals[RESIDUAL_TOTAL]);

        Map<String,Double> transAllocations = geneReadData.getTranscriptAllocations();

        for(int transId = 0; transId < transcriptNames.size(); ++transId)
        {
            double transAllocation = fitAllocations[transId];
            final String trancriptDefn = transcriptNames.get(transId);

            if(transAllocation > 0)
            {
                RE_LOGGER.debug("transcript({}) allocated count({})", trancriptDefn, String.format("%.2f", transAllocation));
            }

            transAllocations.put(trancriptDefn, transAllocation);
        }

        if(mConfig.WriteTransComboData)
        {
            mResultsWriter.writeTransComboCounts(
                    geneReadData, mExpExpressionRates.getCategories(), transComboCounts, fittedCounts);
        }
    }


}
