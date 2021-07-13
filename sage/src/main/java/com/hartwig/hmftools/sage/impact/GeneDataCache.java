package com.hartwig.hmftools.sage.impact;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.impact.ImpactConstants.GENE_UPSTREAM_DISTANCE;

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
        mCurrentGeneIndex = 0;

    }

    public EnsemblDataCache getEnsemblCache() { return mEnsemblDataCache; }

    public boolean loadCache()
    {
        mEnsemblDataCache.setRequiredData(true, false, false, false);
        if(!mEnsemblDataCache.load(false))
            return false;

        mEnsemblDataCache.createGeneIdDataMap();
        return true;
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
        List<EnsemblGeneData> genes = Lists.newArrayList();

        if(mCurrentChromosome == null || !mCurrentChromosome.equals(chromosome))
        {
            mCurrentChromosome = chromosome;
            mCurrentGeneList = mEnsemblDataCache.getChrGeneDataMap().get(chromosome);
            mCurrentGeneIndex = 0;
        }

        if(mCurrentGeneList == null)
        {
            SG_LOGGER.error("invalid chromosome({})", chromosome);
            return genes;
        }

        // setGeneIndex(position);

        for(; mCurrentGeneIndex < mCurrentGeneList.size(); ++mCurrentGeneIndex)
        {
            EnsemblGeneData geneData = mCurrentGeneList.get(mCurrentGeneIndex);

            if(isWithinGeneRange(geneData, position))
            {
                genes.add(geneData);
            }
            else if(position > geneData.GeneEnd)
            {
                continue;
            }
            else if(geneData.GeneStart > position)
            {
                if(mCurrentGeneIndex > 0)
                    --mCurrentGeneIndex;

                break;
            }
        }

        return genes;
    }

    private boolean isWithinGeneRange(final EnsemblGeneData geneData, int position)
    {
        if(geneData.Strand == POS_STRAND)
            return position >= geneData.GeneStart - GENE_UPSTREAM_DISTANCE && position <= geneData.GeneEnd;
        else
            return position >= geneData.GeneStart && position <= geneData.GeneEnd + GENE_UPSTREAM_DISTANCE;
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
