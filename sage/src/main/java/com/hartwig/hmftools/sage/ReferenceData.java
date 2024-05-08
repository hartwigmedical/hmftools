package com.hartwig.hmftools.sage;

import static com.hartwig.hmftools.common.genome.bed.BedFileReader.loadBedFileChrMap;
import static com.hartwig.hmftools.common.hla.HlaCommon.hlaChromosome;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.bed.NamedBedFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.hla.HlaCommon;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotFile;
import com.hartwig.hmftools.sage.common.SimpleVariant;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class ReferenceData
{
    public final ListMultimap<Chromosome,NamedBed> CoveragePanel;
    public final Map<Chromosome,List<BaseRegion>> PanelWithHotspots;
    public final ListMultimap<Chromosome,SimpleVariant> Hotspots;
    public final Map<Chromosome,List<BaseRegion>> HighConfidence;

    public final EnsemblDataCache GeneDataCache;

    public final Map<String,List<TranscriptData>> ChromosomeTranscripts;

    public final IndexedFastaSequenceFile RefGenome;

    private final SageCallConfig mConfig;

    public ReferenceData(final SageCallConfig config, final ConfigBuilder configBuilder)
    {
        mConfig = config;

        CoveragePanel = ArrayListMultimap.create();
        PanelWithHotspots = Maps.newHashMap();
        Hotspots = ArrayListMultimap.create();
        HighConfidence = Maps.newHashMap();
        ChromosomeTranscripts = Maps.newHashMap();

        RefGenome = loadRefGenome(config.Common.RefGenomeFile);

        GeneDataCache = new EnsemblDataCache(configBuilder);
        loadGeneData();

        HlaCommon.populateGeneData(GeneDataCache.getChrGeneDataMap().get(hlaChromosome(config.Common.RefGenVersion)));
    }

    public static IndexedFastaSequenceFile loadRefGenome(final String refGenomeFile)
    {
        try
        {
            return new IndexedFastaSequenceFile(new File(refGenomeFile));
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to load ref genome: {}", e.toString());
            return null;
        }
    }

    private void loadGeneData()
    {
        GeneDataCache.setRequiredData(true, false, false, true);
        GeneDataCache.load(false);

        for(Map.Entry<String,List<GeneData>> entry : GeneDataCache.getChrGeneDataMap().entrySet())
        {
            String chromosome = entry.getKey();

            if(mConfig.Common.SpecificChrRegions.excludeChromosome(chromosome))
                continue;

            List<GeneData> geneDataList = entry.getValue();

            List<TranscriptData> transDataList = Lists.newArrayList();
            ChromosomeTranscripts.put(chromosome, transDataList);

            for(GeneData geneData : geneDataList)
            {
                if(mConfig.Common.SpecificChrRegions.excludeRegion(geneData.GeneStart, geneData.GeneEnd))
                    continue;

                transDataList.add(GeneDataCache.getCanonicalTranscriptData(geneData.GeneId));
            }
        }
    }

    public boolean load()
    {
        if(RefGenome == null)
            return false;

        try
        {
            if(!mConfig.CoverageBed.isEmpty())
            {
                CoveragePanel.putAll(readNamedBedFile(mConfig.CoverageBed));
                SG_LOGGER.info("read {} coverage entries from bed file: {}", CoveragePanel.size(), mConfig.CoverageBed);
            }

            if(!mConfig.PanelBed.isEmpty())
            {
                Map<Chromosome,List<BaseRegion>> panelBed = loadBedFileChrMap(mConfig.PanelBed, true);

                if(panelBed == null)
                    return false;

                PanelWithHotspots.putAll(panelBed);
                SG_LOGGER.info("read {} panel entries from bed file: {}",
                        PanelWithHotspots.values().stream().mapToInt(x -> x.size()).sum(), mConfig.PanelBed);
            }

            loadHotspots();

            if(!mConfig.HighConfidenceBed.isEmpty())
            {
                Map<Chromosome,List<BaseRegion>> hcPanelBed = loadBedFileChrMap(mConfig.HighConfidenceBed);

                if(hcPanelBed == null)
                    return false;

                HighConfidence.putAll(hcPanelBed);
                SG_LOGGER.info("read {} high-confidence entries from bed file: {}",
                        HighConfidence.values().stream().mapToInt(x -> x.size()).sum(), mConfig.HighConfidenceBed);
            }
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to load config");
            return false;
        }

        return true;
    }

    private void loadHotspots() throws IOException
    {
        if(mConfig.Hotspots.isEmpty())
            return;

        ListMultimap<Chromosome,VariantHotspot> hotspotMap = VariantHotspotFile.readFromVCF(mConfig.Hotspots);

        for(VariantHotspot variant : hotspotMap.values())
        {
            Hotspots.put(HumanChromosome.fromString(variant.chromosome()),
                    new SimpleVariant(variant.chromosome(), variant.position(), variant.ref(), variant.alt()));
        }

        SG_LOGGER.info("read {} hotspots from vcf: {}", Hotspots.size(), mConfig.Hotspots);

        // add to panel regions as well if not already covering them
        for(Chromosome chromosome : Hotspots.keySet())
        {
            List<SimpleVariant> hotspots = Hotspots.get(chromosome);
            List<BaseRegion> panelRegions = PanelWithHotspots.get(chromosome);

            if(panelRegions == null)
            {
                panelRegions = Lists.newArrayList();
                PanelWithHotspots.put(chromosome, panelRegions);
            }

            for(SimpleVariant hotspot : hotspots)
            {
                boolean covered = false;
                int index = 0;
                while(index < panelRegions.size())
                {
                    BaseRegion panelRegion = panelRegions.get(index);

                    if(panelRegion.containsPosition(hotspot.position()))
                    {
                        covered = true;
                        break;
                    }
                    else if(hotspot.position() == panelRegion.start() - 1)
                    {
                        covered = true;
                        panelRegion.setStart(hotspot.position());
                        break;
                    }
                    else if(hotspot.position() == panelRegion.end() + 1)
                    {
                        covered = true;
                        panelRegion.setEnd(hotspot.position());

                        if(index < panelRegions.size() - 1)
                        {
                            // check for a merge with the next panel region since BAM slicing does not like adjacent regions
                            BaseRegion nextRegion = panelRegions.get(index + 1);
                            if(nextRegion.start() <= panelRegion.end() + 1)
                            {
                                panelRegion.setEnd(nextRegion.end());
                                panelRegions.remove(index + 1);
                            }
                        }

                        break;
                    }
                    else if(hotspot.position() < panelRegions.get(index).start())
                    {
                        break;
                    }

                    ++index;
                }

                if(!covered)
                {
                    panelRegions.add(index, new BaseRegion(hotspot.position(), hotspot.position()));
                }
            }
        }
    }

    private static ListMultimap<Chromosome,NamedBed> readNamedBedFile(final String panelBed) throws IOException
    {
        final ListMultimap<Chromosome,NamedBed> panel = ArrayListMultimap.create();
        if(!panelBed.isEmpty())
        {
            for(NamedBed bed : NamedBedFile.readBedFile(panelBed))
            {
                if(HumanChromosome.contains(bed.chromosome()))
                {
                    panel.put(HumanChromosome.fromString(bed.chromosome()), bed);
                }
            }
        }

        return panel;
    }
}
