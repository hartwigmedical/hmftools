package com.hartwig.hmftools.sage;

import static com.hartwig.hmftools.common.hla.HlaCommon.hlaChromosome;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.BufferedReader;
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
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotFile;

import org.apache.commons.cli.CommandLine;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class ReferenceData
{
    public final ListMultimap<Chromosome,NamedBed> CoveragePanel;
    public final Map<Chromosome,List<BaseRegion>> PanelWithHotspots;
    public final ListMultimap<Chromosome,VariantHotspot> Hotspots;
    public final Map<Chromosome,List<BaseRegion>> HighConfidence;

    public final EnsemblDataCache GeneDataCache;

    public final Map<String,List<TranscriptData>> ChromosomeTranscripts;

    public final IndexedFastaSequenceFile RefGenome;

    private final SageConfig mConfig;

    public ReferenceData(final SageConfig config, final CommandLine cmd)
    {
        mConfig = config;

        CoveragePanel = ArrayListMultimap.create();
        PanelWithHotspots = Maps.newHashMap();
        Hotspots = ArrayListMultimap.create();
        HighConfidence = Maps.newHashMap();

        RefGenome = loadRefGenome(config.RefGenomeFile);

        GeneDataCache = new EnsemblDataCache(cmd, config.RefGenVersion);
        ChromosomeTranscripts = Maps.newHashMap();
        loadGeneData();

        HlaCommon.populateGeneData(GeneDataCache.getChrGeneDataMap().get(hlaChromosome(config.RefGenVersion)));
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

            if(!mConfig.SpecificChromosomes.isEmpty() && !mConfig.SpecificChromosomes.contains(chromosome))
                continue;

            if(!mConfig.SpecificRegions.isEmpty() && mConfig.SpecificRegions.stream().noneMatch(x -> x.Chromosome.equals(chromosome)))
                continue;

            List<GeneData> geneDataList = entry.getValue();

            List<TranscriptData> transDataList = Lists.newArrayList();
            ChromosomeTranscripts.put(chromosome, transDataList);

            for(GeneData geneData : geneDataList)
            {
                if(!mConfig.SpecificRegions.isEmpty() && mConfig.SpecificRegions.stream()
                        .noneMatch(x -> positionsOverlap(x.start(), x.end(), geneData.GeneStart, geneData.GeneEnd)))
                {
                    continue;
                }

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
                Map<Chromosome,List<BaseRegion>> panelBed = loadBedFile(mConfig.PanelBed);

                if(panelBed == null)
                    return false;

                PanelWithHotspots.putAll(panelBed);
                SG_LOGGER.info("read {} panel entries from bed file: {}",
                        PanelWithHotspots.values().stream().mapToInt(x -> x.size()).sum(), mConfig.PanelBed);
            }

            loadHotspots();

            if(!mConfig.HighConfidenceBed.isEmpty())
            {
                Map<Chromosome,List<BaseRegion>> hcPanelBed = loadBedFile(mConfig.HighConfidenceBed);

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

        Hotspots.putAll(VariantHotspotFile.readFromVCF(mConfig.Hotspots));
        SG_LOGGER.info("read {} hotspots from vcf: {}", Hotspots.size(), mConfig.Hotspots);

        // add to panel regions as well if not already covering them
        for(Chromosome chromosome : Hotspots.keySet())
        {
            List<VariantHotspot> hotspots = Hotspots.get(chromosome);
            List<BaseRegion> panelRegions = PanelWithHotspots.get(chromosome);

            if(panelRegions == null)
            {
                panelRegions = Lists.newArrayList();
                PanelWithHotspots.put(chromosome, panelRegions);
            }

            for(VariantHotspot hotspot : hotspots)
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

    private static Map<Chromosome,List<BaseRegion>> loadBedFile(final String bedFile)
    {
        final Map<Chromosome,List<BaseRegion>> panel = Maps.newHashMap();

        try
        {
            BufferedReader fileReader = createBufferedReader(bedFile);

            final String fileDelim = "\t";
            String line = "";

            List<BaseRegion> chrRegions = null;
            Chromosome currentChr = null;

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(fileDelim, -1);

                Chromosome chromosome = HumanChromosome.fromString(values[0]);
                int posStart = Integer.parseInt(values[1]) + 1; // as per convention
                int posEnd = Integer.parseInt(values[2]);

                if(currentChr != chromosome)
                {
                    currentChr = chromosome;
                    chrRegions = Lists.newArrayList();
                    panel.put(chromosome, chrRegions);
                }

                chrRegions.add(new BaseRegion(posStart, posEnd));
            }
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to load panel BED file({}): {}", bedFile, e.toString());
            return null;
        }

        return panel;
    }
}
