package com.hartwig.hmftools.isofox;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.isofox.IsofoxConfig.GENE_FRAGMENT_BUFFER;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.GeneReadData.markOverlappingGeneRegions;
import static com.hartwig.hmftools.isofox.common.RegionReadData.findUniqueBases;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionsOverlap;
import static com.hartwig.hmftools.isofox.gc.GcRatioCounts.writeReadGcRatioCounts;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.createTranscriptResults;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStr;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.isofox.common.FragmentSizeCalcs;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesData;
import com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesGenerator;
import com.hartwig.hmftools.isofox.exp_rates.ExpectedTransRates;
import com.hartwig.hmftools.isofox.results.GeneResult;
import com.hartwig.hmftools.isofox.results.ResultsWriter;
import com.hartwig.hmftools.isofox.results.TranscriptResult;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;

public class ChromosomeGeneTask implements Callable
{
    private final String mChromosome;
    private final IsofoxConfig mConfig;
    private final SvGeneTranscriptCollection mGeneTransCache;
    private final ResultsWriter mResultsWriter;

    private final GeneBamReader mBamReader;
    private final ExpectedTransRates mExpTransRates;
    private final ExpectedRatesGenerator mExpRatesGenerator;
    private final FragmentSizeCalcs mFragmentSizeCalc;

    private final List<EnsemblGeneData> mGeneDataList;
    private int mCurrentGeneIndex;
    private int mGenesProcessed;

    private final List<GeneResult> mGeneResults;

    private int mCurrentTaskType;
    public static final int CHR_TASK_FRAGMENT_LENGTHS = 0;
    public static final int CHR_TASK_TRANSCRIPT_COUNTS = 1;

    private static final int PERF_TOTAL = 0;
    private static final int PERF_READS = 1;
    private static final int PERF_NOVEL_LOCATIONS = 2;
    private static final int PERF_FIT = 3;
    private final PerformanceCounter[] mPerfCounters;

    public ChromosomeGeneTask(
            final IsofoxConfig config, final String chromosome, final List<EnsemblGeneData> geneDataList,
            final SvGeneTranscriptCollection geneTransCache, final ResultsWriter resultsWriter, final FragmentSizeCalcs fragmentSizeCalc)
    {
        mConfig = config;
        mChromosome = chromosome;
        mGeneTransCache = geneTransCache;
        mResultsWriter = resultsWriter;

        mGeneDataList = geneDataList;

        mCurrentGeneIndex = 0;
        mCurrentTaskType = -1;

        mFragmentSizeCalc = fragmentSizeCalc;

        mBamReader = new GeneBamReader(mConfig, resultsWriter);
        mExpTransRates = mConfig.ApplyExpectedRates ? new ExpectedTransRates(mConfig, resultsWriter) : null;

        mExpRatesGenerator = mConfig.WriteExpectedRates || (mConfig.ApplyExpectedRates && mConfig.UseCalculatedFragmentLengths)
                ? new ExpectedRatesGenerator(mConfig, resultsWriter) : null;

        mGeneResults = Lists.newArrayList();

        mPerfCounters = new PerformanceCounter[PERF_FIT+1];
        mPerfCounters[PERF_TOTAL] = new PerformanceCounter("Total");
        mPerfCounters[PERF_READS] = new PerformanceCounter("ReadCounts");
        mPerfCounters[PERF_NOVEL_LOCATIONS] = new PerformanceCounter("NovelLocations");
        mPerfCounters[PERF_FIT] = new PerformanceCounter("ExpressFit");
    }

    public final GeneBamReader getBamReader() { return mBamReader; }

    public void setTaskType(int taskType) { mCurrentTaskType = taskType; }

    @Override
    public Long call()
    {
        if(mCurrentTaskType == CHR_TASK_TRANSCRIPT_COUNTS)
        {
            assignTranscriptCounts();
            ISF_LOGGER.debug("chromosome({}) transcript counting complete", mChromosome);
        }
        else if(mCurrentTaskType == CHR_TASK_FRAGMENT_LENGTHS)
        {
            calcFragmentLengths();
            ISF_LOGGER.debug("chromosome({}) frag length measurement complete", mChromosome);
        }

        return (long)1; // return value not used
    }

    public void assignTranscriptCounts()
    {
        if(mGeneDataList.size() > 10)
        {
            ISF_LOGGER.info("processing {} genes for chromosome({})", mGeneDataList.size(), mChromosome);
        }

        boolean generateExpRatesOnly = mConfig.WriteExpectedRates && !mConfig.UseCalculatedFragmentLengths && !mConfig.ApplyExpectedRates;

        while(mCurrentGeneIndex < mGeneDataList.size())
        {
            final List<EnsemblGeneData> overlappingGenes = findNextOverlappingGenes();
            final List<GeneReadData> geneReadDataList = createGeneReadData(overlappingGenes);

            for (GeneReadData geneReadData : geneReadDataList)
            {
                mPerfCounters[PERF_TOTAL].start();

                // at the moment it is one or the other
                if(generateExpRatesOnly)
                {
                    generateExpectedTransRates(geneReadData);
                }
                else
                {
                    analyseBamReads(geneReadData);
                }

                mPerfCounters[PERF_TOTAL].stop();

                ISF_LOGGER.debug("chr({}) gene({}) processed({} of {})",
                        mChromosome, geneReadData.name(), mCurrentGeneIndex, mGeneDataList.size());

                ++mGenesProcessed;

                if (mGenesProcessed > 1 && (mGenesProcessed % 100) == 0)
                {
                    ISF_LOGGER.info("chr({}) processed {} of {} genes", mChromosome, mGenesProcessed, mGeneDataList.size());
                }
            }
        }

        writeResults();
    }

    public void calcFragmentLengths()
    {
        mFragmentSizeCalc.calcSampleFragmentSize(mChromosome, mGeneDataList);
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
                ISF_LOGGER.warn("no transcripts found for gene({}:{})", geneData.GeneId, geneData.GeneName);
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
        ISF_LOGGER.info("GENE_OVERLAP: {},{},{},{},{},{}", // chr({}) genes({}) transcripts({}) range({} -> {}),
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
            ISF_LOGGER.warn("invalid gene({}:{}) region({} -> {})", geneData.GeneId, geneData.GeneName, regionStart, regionEnd);
            return;
        }

        GenomeRegion geneRegion = GenomeRegions.create(geneData.Chromosome, regionStart, regionEnd);

        mPerfCounters[PERF_READS].start();
        mBamReader.readBamCounts(geneReadData, geneRegion);
        mPerfCounters[PERF_READS].stop();

        mPerfCounters[PERF_NOVEL_LOCATIONS].start();
        mBamReader.annotateNovelLocations();
        mPerfCounters[PERF_NOVEL_LOCATIONS].stop();

        if(mExpTransRates != null)
        {
            ExpectedRatesData expRatesData = null;

            mPerfCounters[PERF_FIT].start();

            if(mExpRatesGenerator != null)
            {
                generateExpectedTransRates(geneReadData);
                expRatesData = mExpRatesGenerator.getExpectedRatesData();
            }

            mExpTransRates.runTranscriptEstimation(geneReadData, mBamReader.getTransComboData(), expRatesData);

            mPerfCounters[PERF_FIT].stop();
        }

        cacheResults(geneReadData);

        if (mConfig.WriteExonData)
        {
            geneReadData.getTranscripts().forEach(x -> mResultsWriter.writeExonData(geneReadData, x));
        }

        if(mFragmentSizeCalc != null && mConfig.FragmentLengthsByGene)
            mFragmentSizeCalc.writeFragmentLengths(geneData);

        if(mConfig.WriteReadGcRatios)
            writeReadGcRatioCounts(mResultsWriter.getReadGcRatioWriter(), geneData, mBamReader.getGcRatioCounts().getGeneRatioCounts());
    }

    private void cacheResults(final GeneReadData geneReadData)
    {
        final List<TranscriptResult> transResults = Lists.newArrayList();

        if(mConfig.WriteTransData)
        {
            for (final TranscriptData transData : geneReadData.getTranscripts())
            {
                double expRateAllocation = geneReadData.getTranscriptAllocation(transData.TransName);

                final TranscriptResult results =
                        createTranscriptResults(geneReadData, transData, mConfig.ExpRateFragmentLengths, expRateAllocation);

                transResults.add(results);
            }
        }

        GeneResult geneResult = GeneResult.createGeneResults(geneReadData, transResults);

        mGeneResults.add(geneResult);
    }

    private void writeResults()
    {
        for(final GeneResult geneResult : mGeneResults)
        {
            mResultsWriter.writeGeneResult(geneResult);

            geneResult.transcriptResults().forEach(x -> mResultsWriter.writeTranscriptResults(geneResult.geneData(), x));
        }
    }

    public PerformanceCounter[] getPerfCounters()
    {
        return mPerfCounters;
    }


}
