package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.isofox.IsofoxConfig.GENE_FRAGMENT_BUFFER;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RegionReadData.findUniqueBases;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionsOverlap;
import static com.hartwig.hmftools.isofox.gc.GcRatioCounts.writeReadGcRatioCounts;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.createTranscriptResults;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;

import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.isofox.common.FragmentSizeCalcs;
import com.hartwig.hmftools.isofox.common.GeneCollection;
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
    private int mCollectionId;
    private int mCurrentGeneIndex;
    private int mGenesProcessed;

    private final List<GeneResult> mGeneResults;

    private int mCurrentTaskType;
    public static final int CHR_TASK_FRAGMENT_LENGTHS = 0;
    public static final int CHR_TASK_TRANSCRIPT_COUNTS = 1;

    private static final int PERF_TOTAL = 0;
    private static final int PERF_READS = 1;
    private static final int PERF_NOVEL_LOCATIONS = 2;
    public static final int PERF_FIT = 3;
    private static final int PERF_WRITE = 4;
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
        mCollectionId = 0;

        mCurrentGeneIndex = 0;
        mCurrentTaskType = -1;

        mFragmentSizeCalc = fragmentSizeCalc;

        mBamReader = new GeneBamReader(mConfig, resultsWriter);
        mExpTransRates = mConfig.ApplyExpectedRates ? new ExpectedTransRates(mConfig, resultsWriter) : null;

        mExpRatesGenerator = mConfig.writeExpectedRateData() || (mConfig.ApplyExpectedRates && mConfig.UseCalculatedFragmentLengths)
                ? new ExpectedRatesGenerator(mConfig, resultsWriter) : null;

        mGeneResults = Lists.newArrayList();

        mPerfCounters = new PerformanceCounter[PERF_WRITE+1];
        mPerfCounters[PERF_TOTAL] = new PerformanceCounter("Total");
        mPerfCounters[PERF_READS] = new PerformanceCounter("ReadCounts");
        mPerfCounters[PERF_NOVEL_LOCATIONS] = new PerformanceCounter("NovelLocations");
        mPerfCounters[PERF_FIT] = new PerformanceCounter("ExpressFit");

        if(ISF_LOGGER.isInfoEnabled())
            mPerfCounters[PERF_FIT].setSortTimes(true);

        mPerfCounters[PERF_WRITE] = new PerformanceCounter("WriteData");
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

        boolean generateExpRatesOnly = mConfig.generateExpRatesOnly();
        int nextLogCount = 100;

        while(mCurrentGeneIndex < mGeneDataList.size())
        {
            final List<EnsemblGeneData> overlappingGenes = findNextOverlappingGenes();
            final List<GeneReadData> geneReadDataList = createGeneReadData(overlappingGenes);

            GeneCollection geneCollection = new GeneCollection(mCollectionId++, geneReadDataList);

            mPerfCounters[PERF_TOTAL].start();

            // at the moment it is one or the other
            if(generateExpRatesOnly)
            {
                generateExpectedTransRates(geneCollection);
            }
            else
            {
                analyseBamReads(geneCollection);
            }

            mPerfCounters[PERF_TOTAL].stop();

            ISF_LOGGER.debug("chr({}) gene({}) processed({} of {})",
                    mChromosome, geneCollection.geneNames(10), mCurrentGeneIndex, mGeneDataList.size());

            ++mGenesProcessed;

            if (mGenesProcessed >= nextLogCount)
            {
                nextLogCount += 100;
                ISF_LOGGER.info("chr({}) processed {} of {} genes", mChromosome, mGenesProcessed, mGeneDataList.size());
            }
        }

        mPerfCounters[PERF_WRITE].start();
        writeResults();
        mPerfCounters[PERF_WRITE].stop();
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

            if(overlappingGenes.isEmpty()
            || overlappingGenes.stream().anyMatch(x -> positionsOverlap(geneData.GeneStart, geneData.GeneEnd, x.GeneStart, x.GeneEnd)))
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

            geneReadDataList.add(geneReadData);
        }

        return geneReadDataList;
    }

    private void generateExpectedTransRates(final GeneCollection genes)
    {
        mExpRatesGenerator.generateExpectedRates(genes);
    }

    private void analyseBamReads(final GeneCollection geneCollection)
    {
        // cache reference bases for comparison with read bases
        if(mConfig.RefFastaSeqFile != null)
        {
            for (RegionReadData region : geneCollection.getExonRegions())
            {
                final String regionRefBases = mConfig.RefFastaSeqFile.getSubsequenceAt(
                        region.Chromosome, region.PosStart, region.PosEnd).getBaseString();

                region.setRefBases(regionRefBases);
            }

            findUniqueBases(geneCollection.getExonRegions());
        }

        // use a buffer around the gene to pick up reads which span outside its transcripts
        long regionStart = geneCollection.regionBounds()[SE_START] - GENE_FRAGMENT_BUFFER;
        long regionEnd = geneCollection.regionBounds()[SE_END] + GENE_FRAGMENT_BUFFER;

        if(regionStart >= regionEnd)
        {
            ISF_LOGGER.warn("invalid geneCollection(first={} genes={}) region({} -> {})",
                    geneCollection.genes().get(0).name(), regionStart, regionEnd);
            return;
        }

        GenomeRegion geneRegion = GenomeRegions.create(geneCollection.chromosome(), regionStart, regionEnd);

        mPerfCounters[PERF_READS].start();
        mBamReader.readBamCounts(geneCollection, geneRegion);
        mPerfCounters[PERF_READS].stop();

        mPerfCounters[PERF_NOVEL_LOCATIONS].start();
        mBamReader.annotateNovelLocations();
        mPerfCounters[PERF_NOVEL_LOCATIONS].stop();

        if(mExpTransRates != null)
        {
            ExpectedRatesData expRatesData = null;

            if(!mConfig.RunPerfChecks)
            {
                mPerfCounters[PERF_FIT].start();
            }
            else
            {
                int transCount = geneCollection.genes().stream().mapToInt(x -> x.getTranscripts().size()).sum();
                int exonCount = geneCollection.genes().stream().mapToInt(x -> x.getTranscripts().stream().mapToInt(y -> y.exons().size()).sum()).sum();

                final String perfId = String.format("%s genes(%s:%s) trans(%s) exons(%s) range(%d)",
                        geneCollection.chrId(), geneCollection.genes().size(), geneCollection.geneNames(),
                        transCount, exonCount, regionEnd - regionStart);

                mPerfCounters[PERF_FIT].start(perfId);
            }

            if(mExpRatesGenerator != null)
            {
                generateExpectedTransRates(geneCollection);
                expRatesData = mExpRatesGenerator.getExpectedRatesData();
            }

            mExpTransRates.runTranscriptEstimation(geneCollection, mBamReader.getTransComboData(), expRatesData);

            mPerfCounters[PERF_FIT].stop();
        }

        for(GeneReadData geneReadData : geneCollection.genes())
        {
            cacheResults(geneCollection, geneReadData);

            if (mConfig.WriteExonData)
            {
                geneReadData.getTranscripts().forEach(x -> mResultsWriter.writeExonData(geneReadData, x));
            }

            if(mFragmentSizeCalc != null && mConfig.FragmentLengthsByGene)
            {
                mFragmentSizeCalc.writeFragmentLengths(geneReadData.GeneData);
            }

            if(mConfig.WriteReadGcRatios)
            {
                writeReadGcRatioCounts(
                        mResultsWriter.getReadGcRatioWriter(), geneReadData.GeneData, mBamReader.getGcRatioCounts().getGeneRatioCounts());
            }
        }
    }

    private void cacheResults(final GeneCollection geneCollection, final GeneReadData geneReadData)
    {
        final List<TranscriptResult> transResults = Lists.newArrayList();

        if(mConfig.WriteTransData)
        {
            for (final TranscriptData transData : geneReadData.getTranscripts())
            {
                final TranscriptResult results =
                        createTranscriptResults(geneCollection, geneReadData, transData, mConfig.ExpRateFragmentLengths);

                transResults.add(results);
            }
        }

        GeneResult geneResult = GeneResult.createGeneResults(geneCollection, geneReadData, transResults);

        mGeneResults.add(geneResult);
    }

    private void writeResults()
    {
        for(final GeneResult geneResult : mGeneResults)
        {
            mResultsWriter.writeGeneResult(geneResult);

            geneResult.transcriptResults().forEach(x -> mResultsWriter.writeTranscriptResults(geneResult.geneData(), x));
        }

        // written on the fly for now
        // if(mConfig.WriteExpectedRates)
        //    mExpRatesGenerator.writeExpectedRatesData();
    }

    public PerformanceCounter[] getPerfCounters()
    {
        return mPerfCounters;
    }


}
