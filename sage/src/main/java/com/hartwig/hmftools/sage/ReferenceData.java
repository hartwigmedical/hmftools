package com.hartwig.hmftools.sage;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipInputStream;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.SortedSetMultimap;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.bed.NamedBedFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.BEDFileLoader;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotFile;
import com.hartwig.hmftools.sage.config.SageConfig;

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

            if(!mConfig.Chromosomes.isEmpty() && !mConfig.Chromosomes.contains(chromosome))
                continue;

            if(!mConfig.SpecificRegions.isEmpty() && mConfig.SpecificRegions.stream().noneMatch(x -> x.Chromosome.equals(chromosome)))
                continue;

            List<GeneData> geneDataList = entry.getValue();

            if(!mConfig.AppendMode)
            {
                mConfig.Quality.populateGeneData(geneDataList);
            }

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

            loadHotspots();

            /*
            final Map<Chromosome,List<BaseRegion>> panelWithoutHotspots = loadBedFile(mConfig.PanelBed);

            if(!mConfig.PanelBed.isEmpty())
            {
                SG_LOGGER.info("read {} panel entries from bed file: {}", panelWithoutHotspots.size(), mConfig.PanelBed);
            }

            PanelWithHotspots.putAll(panelWithHotspots(panelWithoutHotspots, Hotspots));
            */

            if(!mConfig.PanelBed.isEmpty())
            {
                PanelWithHotspots.putAll(loadBedFile(mConfig.PanelBed));
                SG_LOGGER.info("read {} panel entries from bed file: {}", PanelWithHotspots.size(), mConfig.PanelBed);
            }

            if(!mConfig.HighConfidenceBed.isEmpty())
            {
                HighConfidence.putAll(loadBedFile(mConfig.HighConfidenceBed));
                SG_LOGGER.info("read {} high-confidence entries from bed file: {}", HighConfidence.size(), mConfig.HighConfidenceBed);
            }
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to load config");
            return false;
        }

        // Validate Coverage Bed
        if(mConfig.PanelOnly && !CoveragePanel.isEmpty())
        {
            // TODO: just need to test that the coverage regions are covered by the panel bed
            /*
            if(!GenomeRegionsValidation.isSubset(mPanelWithHotspots.values(), mCoveragePanel.values()))
            {
                throw new IOException("Coverage bed must be a subset of panel bed when running in panel only mode");
            }
            */
        }

        return true;
    }

    private void loadHotspots() throws IOException
    {
        if(!mConfig.Hotspots.isEmpty())
        {
            Hotspots.putAll(VariantHotspotFile.readFromVCF(mConfig.Hotspots));
            SG_LOGGER.info("read {} hotspots from vcf: {}", Hotspots.size(), mConfig.Hotspots);
        }
    }

    /*
    private ListMultimap<Chromosome, ChrBaseRegion> panelWithHotspots(
            final Map<Chromosome,List<BaseRegion>> panelWithoutHotspots, final ListMultimap<Chromosome,VariantHotspot> hotspots)
    {
        final ListMultimap<Chromosome, ChrBaseRegion> result = ArrayListMultimap.create();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            List<ChrBaseRegion> regions = Lists.newArrayList();

            if(panelWithoutHotspots.containsKey(chromosome))
            {
                regions.addAll(panelWithoutHotspots.get(chromosome));
            }

            if(hotspots.containsKey(chromosome))
            {
                hotspots.get(chromosome).forEach(x -> regions.add(new ChrBaseRegion(chromosome.toString(), (int)x.position(), (int)x.end())));
            }

            Collections.sort(regions);

            result.putAll(chromosome, regions);
        }

        return result;
    }
    */

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
            BufferedReader fileReader = bedFile.endsWith(".gz") ?
                    new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(bedFile)))) :
                    new BufferedReader(new FileReader(bedFile));

            final String fileDelim = "\t";
            String line = "";

            List<BaseRegion> chrRegions = null;
            Chromosome currentChr = null;

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(fileDelim, -1);

                Chromosome chromosome = HumanChromosome.fromString(values[0]);
                int posStart = Integer.parseInt(values[1]);
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
        }

        return panel;
    }
}
