package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.genome.genepanel.HmfTranscriptRegionFile.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.PaveConstants.GENE_UPSTREAM_DISTANCE;
import static com.hartwig.hmftools.pave.PaveUtils.withinTransRange;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.GeneMappingData;
import com.hartwig.hmftools.common.ensemblcache.GeneNameMapping;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class GeneDataCache
{
    private final EnsemblDataCache mEnsemblDataCache;

    private final List<String> mDriverGenes;
    private final Map<String,List<String>> mOtherReportableTranscripts;

    private final Set<String> mMissingTranscripts; // no longer in Ensembl
    private final Map<String,TranscriptData> mTranscriptDataMap; // transcripts keyed by trans name, built on demand

    private final GeneNameMapping mGeneNameMapping; // only relevant when converting from or checking SnpEff names

    private final boolean mUseIndexing;
    private String mCurrentChromosome;
    private List<GeneData> mCurrentChromosomeGenes;
    private List<GeneData> mCurrentGenes; // in the current vacinity
    private int mCurrentPosStrandGeneIndex;
    private int mCurrentNegStrandGeneIndex;

    public GeneDataCache(
            final String ensemblDir, final RefGenomeVersion refGenVersion, final String driverGeneFile,
            boolean requireMapping, boolean useIndexing)
    {
        mEnsemblDataCache = new EnsemblDataCache(ensemblDir, refGenVersion);

        mDriverGenes = Lists.newArrayList();
        mOtherReportableTranscripts = Maps.newHashMap();
        loadDriverGenes(driverGeneFile);

        mTranscriptDataMap = Maps.newHashMap();
        mMissingTranscripts = Sets.newHashSet();

        mUseIndexing = useIndexing;
        mCurrentChromosome = null;
        mCurrentChromosomeGenes = null;
        mCurrentPosStrandGeneIndex = 0;
        mCurrentNegStrandGeneIndex = 0;
        mCurrentGenes = Lists.newArrayList();

        mGeneNameMapping = requireMapping ? new GeneNameMapping() : null;
    }

    public EnsemblDataCache getEnsemblCache() { return mEnsemblDataCache; }
    public List<String> getDriverPanelGenes() { return mDriverGenes; }
    public Map<String,List<String>> getOtherReportableTranscripts() { return mOtherReportableTranscripts; }

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
                    mOtherReportableTranscripts.put(driverGene.gene(), driverGene.additionalReportedTranscripts());
            }

            if(!mOtherReportableTranscripts.isEmpty())
            {
                PV_LOGGER.info("loaded {} driver alternative transcripts from {} genes",
                        mOtherReportableTranscripts.values().stream().mapToInt(x -> x.size()).sum(), mOtherReportableTranscripts.size());
            }
        }
        catch (IOException e)
        {
            PV_LOGGER.error("failed to load driver gene panel file({}): {}", driverGeneFile, e.toString());
        }
    }

    public List<TranscriptData> findTranscripts(final String geneId, int startPosition, int endPosition)
    {
        List<TranscriptData> transDataList = mEnsemblDataCache.getTranscripts(geneId);

        if(transDataList == null)
            return Lists.newArrayList();

        return transDataList.stream()
                .filter(x -> withinTransRange(x, startPosition, endPosition))
                .collect(Collectors.toList());
    }

    public TranscriptData getTranscriptData(final String geneId, final String transName)
    {
        if(mMissingTranscripts.contains(transName))
            return null;

        TranscriptData transData = mTranscriptDataMap.get(transName);

        if(transData != null)
            return transData;

        transData = mEnsemblDataCache.getTranscriptData(geneId, transName);

        if(transData == null)
        {
            mMissingTranscripts.add(transName);
            return null;
        }

        mTranscriptDataMap.put(transName, transData);
        return transData;
    }

    public List<GeneData> findGenes(final String chromosome, int startPosition, int endPosition)
    {
        if(!mUseIndexing)
        {
            List<GeneData> genes = Lists.newArrayList();

            List<GeneData> geneDataList = mEnsemblDataCache.getChrGeneDataMap().get(chromosome);

            if(geneDataList == null)
                return genes;

            for(GeneData geneData : geneDataList)
            {
                if(isWithinGeneRange(geneData, startPosition, endPosition))
                    genes.add(geneData);
                else if(geneData.GeneStart > endPosition + GENE_UPSTREAM_DISTANCE)
                    break;
            }

            return genes;
        }

        if(mCurrentChromosome == null || !mCurrentChromosome.equals(chromosome))
        {
            mCurrentChromosome = chromosome;
            mCurrentChromosomeGenes = mEnsemblDataCache.getChrGeneDataMap().get(chromosome);
            mCurrentPosStrandGeneIndex = 0;
            mCurrentNegStrandGeneIndex = 0;
            mCurrentGenes.clear();
        }

        if(mCurrentChromosomeGenes == null)
        {
            // PV_LOGGER.error("invalid chromosome({})", chromosome); // just MT
            return Lists.newArrayList();
        }

        // purge any gene where the position is now past its end
        int index = 0;
        while(index < mCurrentGenes.size())
        {
            GeneData geneData = mCurrentGenes.get(index);

            boolean isPastGene = (geneData.forwardStrand() && geneData.GeneEnd < endPosition)
                    || (geneData.reverseStrand() && geneData.GeneEnd + GENE_UPSTREAM_DISTANCE < endPosition);

            if(isPastGene)
                mCurrentGenes.remove(index);
            else
                ++index;
        }

        // otherwise search forward and add any additional genes which overlap the position
        // need to maintain 2 indices since the range checks can result in different current positions
        for(int i = 0; i <= 1; ++i)
        {
            int currentStrand = (i == 0) ? POS_STRAND : NEG_STRAND;
            int geneIndex = (i == 0) ? mCurrentPosStrandGeneIndex : mCurrentNegStrandGeneIndex;

            for(; geneIndex < mCurrentChromosomeGenes.size(); ++geneIndex)
            {
                GeneData geneData = mCurrentChromosomeGenes.get(geneIndex);

                if(geneData.Strand != currentStrand)
                    continue;

                if(isWithinGeneRange(geneData, startPosition, endPosition))
                {
                    if(mCurrentGenes.contains(geneData))
                    {
                        PV_LOGGER.error("adding current gene({}:{}) index({}) twice", geneData.GeneId, geneData.GeneName, geneIndex);
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

            if(i == 0)
                mCurrentPosStrandGeneIndex = geneIndex;
            else
                mCurrentNegStrandGeneIndex = geneIndex;
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

    public GeneData findSnpEffGeneData(final String geneId, final String geneName)
    {
        if(geneId.isEmpty())
            return null;

        GeneData geneData = mEnsemblDataCache.getGeneDataById(geneId);

        if(geneData != null)
            return geneData;

        GeneMappingData mappingData = mGeneNameMapping.getMappingDataByOld(geneName);

        if(mappingData == null)
            return null;

        return mEnsemblDataCache.getGeneDataById(mappingData.GeneId);
    }

    public String getSnpEffGeneName(final String geneName)
    {
        return mGeneNameMapping.isUnchanged(geneName) ? geneName : mGeneNameMapping.getOldName(geneName);
    }

    public String getGeneNameFromSnpEff(final String snpEffGeneName)
    {
        return mGeneNameMapping.getNewName(snpEffGeneName);
    }
}
