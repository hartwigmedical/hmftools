package com.hartwig.hmftools.sage.impact;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.impact.ImpactConstants.GENE_UPSTREAM_DISTANCE;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.compress.utils.Lists;

public class GeneDataCache
{
    private final EnsemblDataCache mEnsemblDataCache;

    private final List<String> mDriverGenes;

    private String mCurrentChromosome;
    private List<GeneData> mCurrentChromosomeGenes;
    private List<GeneData> mCurrentGenes; // in the current vacinity
    private int mCurrentGeneIndex;

    public GeneDataCache(final String ensemblDir, final RefGenomeVersion refGenomeVersion, final String driverGenePanel)
    {
        mEnsemblDataCache = new EnsemblDataCache(ensemblDir, refGenomeVersion);

        mDriverGenes = Lists.newArrayList();
        loadDriverGenes(driverGenePanel);

        mCurrentChromosome = null;
        mCurrentChromosomeGenes = null;
        mCurrentGeneIndex = 0;
        mCurrentGenes = Lists.newArrayList();
    }

    public EnsemblDataCache getEnsemblCache() { return mEnsemblDataCache; }
    public List<String> getDriverPanelGenes() { return mDriverGenes; }

    public boolean loadCache()
    {
        mEnsemblDataCache.setRequiredData(true, false, false, false);
        if(!mEnsemblDataCache.load(false))
            return false;

        mEnsemblDataCache.createGeneIdDataMap();
        return true;
    }

    private void loadDriverGenes(final String driverGeneFile)
    {
        if(driverGeneFile == null || driverGeneFile.isEmpty())
            return;

        try
        {
            List<DriverGene> driverGenes = DriverGeneFile.read(driverGeneFile);
            driverGenes.forEach(x -> mDriverGenes.add(x.gene()));
        }
        catch (IOException e)
        {
            SG_LOGGER.error("failed to load driver gene panel file({}): {}", driverGeneFile, e.toString());
        }
    }

    public List<TranscriptData> findTranscripts(final String geneId, int position)
    {
        List<TranscriptData> transDataList = mEnsemblDataCache.getTranscripts(geneId);

        return transDataList.stream()
                .filter(x -> x.TransStart <= position && position <= x.TransEnd)
                .filter(x -> x.CodingStart != null)
                .collect(Collectors.toList());
    }

    public List<GeneData> findGenes(final String chromosome, int position)
    {
        if(mCurrentChromosome == null || !mCurrentChromosome.equals(chromosome))
        {
            mCurrentChromosome = chromosome;
            mCurrentChromosomeGenes = mEnsemblDataCache.getChrGeneDataMap().get(chromosome);
            mCurrentGeneIndex = 0;
            mCurrentGenes.clear();
        }

        if(mCurrentChromosomeGenes == null)
        {
            SG_LOGGER.error("invalid chromosome({})", chromosome);
            return Lists.newArrayList();
        }

        // purge any gene where the position is now past its end
        int index = 0;
        while(index < mCurrentGenes.size())
        {
            if(mCurrentGenes.get(index).GeneEnd < position)
                mCurrentGenes.remove(index);
            else
                ++index;
        }

        for(; mCurrentGeneIndex < mCurrentChromosomeGenes.size(); ++mCurrentGeneIndex)
        {
            GeneData geneData = mCurrentChromosomeGenes.get(mCurrentGeneIndex);

            if(isWithinGeneRange(geneData, position))
            {
                if(mCurrentGenes.contains(geneData))
                {
                    SG_LOGGER.error("adding current gene({}:{}) index({}) twice", geneData.GeneId, geneData.GeneName, mCurrentGeneIndex);
                    break;
                }

                mCurrentGenes.add(geneData);
            }
            else if(position > geneData.GeneEnd)
            {
                continue;
            }
            else if(geneData.GeneStart > position)
            {
                break;
            }
        }

        return mCurrentGenes;
    }

    private boolean isWithinGeneRange(final GeneData geneData, int position)
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
