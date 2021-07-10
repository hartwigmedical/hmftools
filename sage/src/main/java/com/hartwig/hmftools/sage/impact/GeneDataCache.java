package com.hartwig.hmftools.sage.impact;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.compress.utils.Lists;

public class GeneDataCache
{
    private final EnsemblDataCache mEnsemblDataCache;

    private String mCurrentChromosome;
    private List<EnsemblGeneData> mCurrentGeneList;
    private int mCurrentGeneIndex;

    public GeneDataCache(final String ensemblDir, final RefGenomeVersion refGenomeVersion)
    {
        mEnsemblDataCache = new EnsemblDataCache(ensemblDir, refGenomeVersion);

        mCurrentChromosome = null;
        mCurrentGeneList = null;
        mCurrentGeneIndex = -1;

    }

    public EnsemblDataCache getEnsemblCache() { return mEnsemblDataCache; }

    public void loadCache()
    {
        mEnsemblDataCache.setRequiredData(true, false, false, false);
        mEnsemblDataCache.load(false);
    }

    public List<TranscriptData> findTranscripts(final String geneId, int position)
    {
        List<TranscriptData> transDataList = mEnsemblDataCache.getTranscripts(geneId);

        return transDataList.stream()
                .filter(x -> x.TransStart <= position && position <= x.TransEnd)
                .filter(x -> x.CodingStart != null)
                .collect(Collectors.toList());
    }

    public List<EnsemblGeneData> findGenes(final String chromosome, int position)
    {
        if(mCurrentChromosome == null || !mCurrentChromosome.equals(chromosome))
        {
            mCurrentChromosome = chromosome;
            mCurrentGeneList = mEnsemblDataCache.getChrGeneDataMap().get(chromosome);
            mCurrentGeneIndex = 0;
        }

        // setGeneIndex(position);

        List<EnsemblGeneData> genes = Lists.newArrayList();

        for(; mCurrentGeneIndex < mCurrentGeneList.size(); ++mCurrentGeneIndex)
        {
            EnsemblGeneData geneData = mCurrentGeneList.get(mCurrentGeneIndex);

            if(position > geneData.GeneEnd)
            {
                ++mCurrentGeneIndex;
                continue;
            }

            if(geneData.GeneStart <= position && geneData.GeneEnd >= position)
            {
                genes.add(geneData);
            }
            else if(geneData.GeneStart > position)
            {
                --mCurrentGeneIndex;
                break;
            }
        }

        return genes;
    }

    /*
    private void setGeneIndex(int position)
    {
        // move the current index to cover the position or be the last gene before it
        for(; mCurrentGeneIndex < mCurrentGeneList.size(); ++mCurrentGeneIndex)
        {
            EnsemblGeneData geneData = mCurrentGeneList.get(mCurrentGeneIndex);

            if(geneData.GeneStart > position)
            {
                --mCurrentGeneIndex;
                return;
            }

            if(position > geneData.GeneEnd)
            {
                ++mCurrentGeneIndex;
            }
        }
    }
    */
}
