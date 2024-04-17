package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.PaveConstants.GENE_UPSTREAM_DISTANCE;
import static com.hartwig.hmftools.pave.impact.PaveUtils.withinTransRange;

import java.io.IOException;
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
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.jetbrains.annotations.Nullable;

public class GeneDataCache
{
    private final EnsemblDataCache mEnsemblDataCache;

    private final String mDriverGeneFile;
    private final List<DriverGene> mDriverGenes;
    private final Set<String> mDriverGeneNames;
    private final Map<String,List<String>> mOtherReportableTranscripts;

    public GeneDataCache(final String ensemblDir, final RefGenomeVersion refGenVersion, final String driverGeneFile)
    {
        mEnsemblDataCache = new EnsemblDataCache(ensemblDir, refGenVersion);

        mDriverGeneFile = driverGeneFile;
        mDriverGenes = Lists.newArrayList();
        mDriverGeneNames = Sets.newHashSet();
        mOtherReportableTranscripts = Maps.newHashMap();
    }

    public EnsemblDataCache getEnsemblCache() { return mEnsemblDataCache; }
    public boolean isDriverPanelGene(final String geneName) { return mDriverGeneNames.contains(geneName); }
    public List<DriverGene> getDriverPanel() { return mDriverGenes; }
    public Map<String,List<String>> getOtherReportableTranscripts() { return mOtherReportableTranscripts; }

    public boolean loadCache(boolean canonicalOnly, boolean onlyDriverGenes)
    {
        if(!loadDriverGenes(mDriverGeneFile))
            return false;

        mEnsemblDataCache.setRequiredData(true, false, false, canonicalOnly);
        mEnsemblDataCache.setRequireNonEnsemblTranscripts();

        if(onlyDriverGenes)
        {
            if(!mEnsemblDataCache.load(true))
            {
                PV_LOGGER.error("invalid Ensembl data cache");
                return false;
            }

            for(String geneName : mDriverGeneNames)
            {
                GeneData geneData = mEnsemblDataCache.getGeneDataByName(geneName);

                if(geneData == null)
                {
                    PV_LOGGER.error("gene({}) missing from Ensembl cache", geneName);
                    return false;
                }
            }

            List<String> driverGeneIds = mDriverGeneNames.stream()
                    .map(x -> mEnsemblDataCache.getGeneDataByName(x))
                    .filter(x -> x != null)
                    .map(x -> x.GeneId)
                    .collect(Collectors.toList());

            mEnsemblDataCache.loadTranscriptData(driverGeneIds);
        }
        else
        {
            if(!mEnsemblDataCache.load(false))
            {
                PV_LOGGER.error("invalid Ensembl data cache");
                return false;
            }
        }

        mEnsemblDataCache.createGeneIdDataMap();

        return true;
    }

    private boolean loadDriverGenes(final String driverGeneFile)
    {
        if(driverGeneFile == null || driverGeneFile.isEmpty())
            return true;

        try
        {
            mDriverGenes.addAll(DriverGeneFile.read(driverGeneFile));

            for(DriverGene driverGene : mDriverGenes)
            {
                mDriverGeneNames.add(driverGene.gene());

                if(!driverGene.additionalReportedTranscripts().isEmpty())
                    mOtherReportableTranscripts.put(driverGene.gene(), driverGene.additionalReportedTranscripts());
            }

            if(!mOtherReportableTranscripts.isEmpty())
            {
                PV_LOGGER.debug("loaded {} driver alternative transcripts from {} genes",
                        mOtherReportableTranscripts.values().stream().mapToInt(x -> x.size()).sum(), mOtherReportableTranscripts.size());
            }
        }
        catch (IOException e)
        {
            PV_LOGGER.error("failed to load driver gene panel file({}): {}", driverGeneFile, e.toString());
            return false;
        }

        return true;
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

    public GeneCacheIndexing createIndexing(final String chromosome)
    {
        return new GeneCacheIndexing(mEnsemblDataCache.getChrGeneDataMap().get(chromosome));
    }

    public List<GeneData> findGenes(
            final String chromosome, int startPosition, int endPosition, @Nullable final GeneCacheIndexing cacheIndexing)
    {
        if(cacheIndexing == null)
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

        if(cacheIndexing.ChromosomeGenes == null)
        {
            // PV_LOGGER.error("invalid chromosome({})", chromosome); // just MT
            return Lists.newArrayList();
        }

        // purge any gene where the position is now past its end
        int index = 0;
        while(index < cacheIndexing.CurrentGenes.size())
        {
            GeneData geneData = cacheIndexing.CurrentGenes.get(index);

            boolean isPastGene = (geneData.forwardStrand() && geneData.GeneEnd < endPosition)
                    || (geneData.reverseStrand() && geneData.GeneEnd + GENE_UPSTREAM_DISTANCE < endPosition);

            if(isPastGene)
                cacheIndexing.CurrentGenes.remove(index);
            else
                ++index;
        }

        // otherwise search forward and add any additional genes which overlap the position
        // need to maintain 2 indices since the range checks can result in different current positions
        for(int i = 0; i <= 1; ++i)
        {
            int currentStrand = (i == 0) ? POS_STRAND : NEG_STRAND;
            int geneIndex = (i == 0) ? cacheIndexing.CurrentPosStrandGeneIndex : cacheIndexing.CurrentNegStrandGeneIndex;

            for(; geneIndex < cacheIndexing.ChromosomeGenes.size(); ++geneIndex)
            {
                GeneData geneData = cacheIndexing.ChromosomeGenes.get(geneIndex);

                if(geneData.Strand != currentStrand)
                    continue;

                if(isWithinGeneRange(geneData, startPosition, endPosition))
                {
                    if(cacheIndexing.CurrentGenes.contains(geneData))
                    {
                        PV_LOGGER.error("adding current gene({}:{}) index({}) twice", geneData.GeneId, geneData.GeneName, geneIndex);
                        break;
                    }

                    cacheIndexing.CurrentGenes.add(geneData);
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
                cacheIndexing.CurrentPosStrandGeneIndex = geneIndex;
            else
                cacheIndexing.CurrentNegStrandGeneIndex = geneIndex;
        }

        return cacheIndexing.CurrentGenes;
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
