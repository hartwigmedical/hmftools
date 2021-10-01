package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.genome.genepanel.HmfTranscriptRegionFile.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.pave.PaveConfig.VI_LOGGER;
import static com.hartwig.hmftools.pave.PaveConstants.GENE_UPSTREAM_DISTANCE;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
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
    private final Map<String,List<String>> mOtherReportableTranscripts;

    private String mCurrentChromosome;
    private List<GeneData> mCurrentChromosomeGenes;
    private List<GeneData> mCurrentGenes; // in the current vacinity
    private int mCurrentGeneIndex;

    public GeneDataCache(final String ensemblDir, final RefGenomeVersion refGenomeVersion, final String driverGenePanel)
    {
        mEnsemblDataCache = new EnsemblDataCache(ensemblDir, refGenomeVersion);

        mDriverGenes = Lists.newArrayList();
        mOtherReportableTranscripts = Maps.newHashMap();
        loadDriverGenes(driverGenePanel);

        mCurrentChromosome = null;
        mCurrentChromosomeGenes = null;
        mCurrentGeneIndex = 0;
        mCurrentGenes = Lists.newArrayList();
    }

    public EnsemblDataCache getEnsemblCache() { return mEnsemblDataCache; }
    public List<String> getDriverPanelGenes() { return mDriverGenes; }
    public Map<String,List<String>> getOtherReportableTranscripts() { return mOtherReportableTranscripts; }

    public boolean loadCache() { return loadCache(false, false); }

    public boolean loadCache(boolean canonicalOnly, boolean onlyDriverGenes)
    {
        mEnsemblDataCache.setRequiredData(true, false, false, canonicalOnly);

        if(onlyDriverGenes)
        {
            if(!mEnsemblDataCache.load(true))
                return false;

            List<String> driverGeneIds = mDriverGenes.stream()
                    .map(x -> mEnsemblDataCache.getGeneDataByName(x).GeneId).collect(Collectors.toList());

            mEnsemblDataCache.loadTranscriptData(driverGeneIds);
        }
        else
        {
            if(!mEnsemblDataCache.load(false))
                return false;
        }

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

            for(DriverGene driverGene : driverGenes)
            {
                driverGenes.forEach(x -> mDriverGenes.add(x.gene()));

                if(!driverGene.additionalReportedTranscripts().isEmpty())
                {
                    mOtherReportableTranscripts.put(
                            driverGene.gene(),
                            Arrays.stream(driverGene.additionalReportedTranscripts().split(ITEM_DELIM)).collect(Collectors.toList()));
                }
            }

            if(!mOtherReportableTranscripts.isEmpty())
            {
                VI_LOGGER.info("loaded {} driver alternative transcripts from {} genes",
                        mOtherReportableTranscripts.values().stream().mapToInt(x -> x.size()).sum(), mOtherReportableTranscripts.size());
            }
        }
        catch (IOException e)
        {
            VI_LOGGER.error("failed to load driver gene panel file({}): {}", driverGeneFile, e.toString());
        }
    }

    public List<TranscriptData> findTranscripts(final String geneId, int startPosition, int endPosition)
    {
        List<TranscriptData> transDataList = mEnsemblDataCache.getTranscripts(geneId);

        return transDataList.stream()
                .filter(x -> positionsOverlap(startPosition, endPosition, x.TransStart, x.TransEnd))
                .filter(x -> x.CodingStart != null)
                .collect(Collectors.toList());
    }

    public List<GeneData> findGenes(final String chromosome, int startPosition, int endPosition)
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
            VI_LOGGER.error("invalid chromosome({})", chromosome);
            return Lists.newArrayList();
        }

        // purge any gene where the position is now past its end
        int index = 0;
        while(index < mCurrentGenes.size())
        {
            if(mCurrentGenes.get(index).GeneEnd < endPosition)
                mCurrentGenes.remove(index);
            else
                ++index;
        }

        // otherwise search forward and add any additional genes which overlap the position
        for(; mCurrentGeneIndex < mCurrentChromosomeGenes.size(); ++mCurrentGeneIndex)
        {
            GeneData geneData = mCurrentChromosomeGenes.get(mCurrentGeneIndex);

            if(isWithinGeneRange(geneData, startPosition, endPosition))
            {
                if(mCurrentGenes.contains(geneData))
                {
                    VI_LOGGER.error("adding current gene({}:{}) index({}) twice", geneData.GeneId, geneData.GeneName, mCurrentGeneIndex);
                    break;
                }

                mCurrentGenes.add(geneData);
            }
            else if(startPosition > geneData.GeneEnd)
            {
                continue;
            }
            else if(geneData.GeneStart > startPosition)
            {
                break;
            }
        }

        return mCurrentGenes;
    }

    private boolean isWithinGeneRange(final GeneData geneData, int startPosition, int endPosition)
    {
        int geneRangeStart = geneData.GeneStart;
        int geneRangeEnd = geneData.GeneEnd;

        if(geneData.Strand == POS_STRAND)
            geneRangeStart -= GENE_UPSTREAM_DISTANCE;
        else
            geneRangeEnd += GENE_UPSTREAM_DISTANCE;

        return positionsOverlap(startPosition, endPosition, geneRangeStart, geneRangeEnd);
    }
}
