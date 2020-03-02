package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStr;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.markOverlappingGeneRegions;
import static com.hartwig.hmftools.svtools.rna_expression.RegionReadData.findUniqueBases;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.GENE_FRAGMENT_BUFFER;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.RE_LOGGER;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpUtils.positionsOverlap;
import static com.hartwig.hmftools.svtools.rna_expression.TranscriptModel.calculateTranscriptResults;

import java.util.List;
import java.util.concurrent.Callable;
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

    private final GeneBamReader mBamReader;
    private final ExpectedTransRates mExpTransRates;
    private final ExpectedRatesGenerator mExpRatesGenerator;
    private final FragmentSizeCalcs mFragmentSizeCalc;

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

        mFragmentSizeCalc = mConfig.WriteFragmentLengths || mConfig.UseCalculatedFragmentLengths
                ? new FragmentSizeCalcs(config, null, resultsWriter.getFragmentLengthWriter()) : null;

        mBamReader = new GeneBamReader(mConfig, resultsWriter, mFragmentSizeCalc);
        mExpTransRates = mConfig.ApplyExpectedRates ? new ExpectedTransRates(mConfig, resultsWriter) : null;

        mExpRatesGenerator = mConfig.WriteExpectedRates || (mConfig.ApplyExpectedRates && mConfig.UseCalculatedFragmentLengths)
                ? new ExpectedRatesGenerator(mConfig, resultsWriter) : null;

        mPerfCounter = new PerformanceCounter(String.format("chr(%s) genes(%d)", mChromosome, mGeneDataList.size()));
    }

    public final GeneBamReader getBamReader() { return mBamReader; }

    @Override
    public Long call()
    {
        analyseGenes();
        return (long)1; // return value not used
    }

    public void analyseGenes()
    {
        if(mGeneDataList.size() > 10)
        {
            RE_LOGGER.info("processing {} genes for chromosome({})", mGeneDataList.size(), mChromosome);
        }

        boolean generateExpRatesOnly = mConfig.WriteExpectedRates && !mConfig.UseCalculatedFragmentLengths && !mConfig.ApplyExpectedRates;

        while(mCurrentGeneIndex < mGeneDataList.size())
        {
            final List<EnsemblGeneData> overlappingGenes = findNextOverlappingGenes();
            final List<GeneReadData> geneReadDataList = createGeneReadData(overlappingGenes);

            for (GeneReadData geneReadData : geneReadDataList)
            {
                mPerfCounter.start();

                RE_LOGGER.debug("chr({}) gene({}) processed({} of {})",
                        mChromosome, geneReadData.name(), mCurrentGeneIndex, mGeneDataList.size());

                // at the moment it is one or the other
                if(generateExpRatesOnly)
                {
                    generateExpectedTransRates(geneReadData);
                }
                else
                {
                    analyseBamReads(geneReadData);
                }

                mPerfCounter.stop();

                ++mGenesProcessed;

                if (mGenesProcessed > 1 && (mGenesProcessed % 100) == 0)
                {
                    RE_LOGGER.info("chr({}) processed {} genes", mChromosome, mGenesProcessed);
                }
            }
        }
    }

    public PerformanceCounter getPerfStats()
    {
        return mPerfCounter;
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

        markOverlappingGeneRegions(geneReadDataList, false);

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
            geneNamesStr = appendStr(geneNamesStr, geneReadData.GeneData.GeneId, ';');
            transcriptCount += geneReadData.getTranscripts().size();
            maxRange =  max(maxRange, geneReadData.GeneData.GeneEnd);
            minRange =  minRange < 0 ? geneReadData.GeneData.GeneStart : min(minRange, geneReadData.GeneData.GeneStart);
        }

        // Time,Chromosome,GeneCount,TranscriptCount,RangeStart,RangeEnd,GeneNames
        RE_LOGGER.info("GENE_OVERLAP: {},{},{},{},{},{}", // chr({}) genes({}) transcripts({}) range({} -> {}),
                mChromosome, overlappingGenes.size(), transcriptCount, minRange, maxRange, geneNamesStr);
    }

    private void generateExpectedTransRates(final GeneReadData geneReadData)
    {
        mExpRatesGenerator.generateExpectedRates(geneReadData);
    }

    private void analyseBamReads(final GeneReadData geneReadData)
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

        mBamReader.readBamCounts(geneReadData, geneRegion);

        if(mExpTransRates != null)
        {
            ExpectedRatesData expRatesData = null;

            if(mExpRatesGenerator != null)
            {
                if (mConfig.UseCalculatedFragmentLengths)
                    mFragmentSizeCalc.setConfigLengthDistribution();

                generateExpectedTransRates(geneReadData);
                expRatesData = mExpRatesGenerator.getExpectedRatesData();
            }

            mExpTransRates.runTranscriptEstimation(geneReadData, mBamReader.getTransComboData(), expRatesData);
        }

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

        if(mFragmentSizeCalc != null)
            mFragmentSizeCalc.writeFragmentLengths(geneData);
    }

}
