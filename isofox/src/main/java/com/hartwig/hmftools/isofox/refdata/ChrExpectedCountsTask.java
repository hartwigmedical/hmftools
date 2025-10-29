package com.hartwig.hmftools.isofox.refdata;

import static com.hartwig.hmftools.isofox.ChromosomeTaskExecutor.findNextOverlappingGenes;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.GeneReadData.createGeneReadData;

import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;

public class ChrExpectedCountsTask implements Callable<Void>
{
    private final RefDataConfig mConfig;
    private final EnsemblDataCache mGeneTransCache;

    private final ExpectedCountsGenerator mExpRatesGenerator;

    private String mChromosome;
    private final List<GeneData> mGeneDataList;

    private int mCurrentGeneIndex;
    private int mCollectionId;
    private int mGenesProcessed;

    public ChrExpectedCountsTask(final RefDataConfig config, final EnsemblDataCache ensemblDataCache, final RefDataWriter writer)
    {
        mConfig = config;
        mGeneTransCache = ensemblDataCache;

        mGeneDataList = Lists.newArrayList();
        mChromosome = "";
        mCurrentGeneIndex = 0;

        mExpRatesGenerator = new ExpectedCountsGenerator(mConfig, writer);
    }

    public void initialise(final String chromosome, final List<GeneData> geneDataList)
    {
        mChromosome = chromosome;
        mGeneDataList.clear();
        mGeneDataList.addAll(geneDataList);
        mCurrentGeneIndex = 0;
    }

    @Override
    public Void call()
    {
        generateExpectedCounts();
        return null;
    }

    private void generateExpectedCounts()
    {
        if(mGeneDataList.size() > 10)
        {
            ISF_LOGGER.info("processing {} genes for chromosome({})", mGeneDataList.size(), mChromosome);
        }

        mCurrentGeneIndex = 0;
        final List<GeneData> overlappingGenes = Lists.newArrayList();
        int nextLogCount = 100;

        while(mCurrentGeneIndex < mGeneDataList.size())
        {
            mCurrentGeneIndex = findNextOverlappingGenes(mGeneDataList, mCurrentGeneIndex, overlappingGenes);
            final List<GeneReadData> geneReadDataList = createGeneReadData(overlappingGenes, mGeneTransCache);

            GeneCollection geneCollection = new GeneCollection(mCollectionId++, geneReadDataList);

            for(GeneReadData geneReadData : geneReadDataList)
            {
                if(mConfig.EnrichedGeneIds.contains(geneReadData.GeneData.GeneId))
                {
                    geneCollection.setEnrichedTranscripts(mGeneTransCache.getTranscripts(geneReadData.GeneData.GeneId));
                }
            }

            mExpRatesGenerator.generateExpectedRates(geneCollection);

            ISF_LOGGER.trace("chr({}) gene({}) processed({} of {})",
                    mChromosome, geneCollection.geneNames(10), mCurrentGeneIndex, mGeneDataList.size());

            ++mGenesProcessed;

            if(mGenesProcessed >= nextLogCount)
            {
                nextLogCount += 100;
                ISF_LOGGER.debug("chr({}) processed {} of {} genes", mChromosome, mGenesProcessed, mGeneDataList.size());
            }
        }

        if(nextLogCount > 100)
            ISF_LOGGER.info("chromosome({}) expected transcript counts complete", mChromosome);
    }
}
