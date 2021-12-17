package com.hartwig.hmftools.sage;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.SortedSetMultimap;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.bed.NamedBedFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.region.BEDFileLoader;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotFile;
import com.hartwig.hmftools.sage.config.SageConfig;

import org.apache.commons.compress.utils.Lists;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class ReferenceData
{
    public final ListMultimap<Chromosome,NamedBed> CoveragePanel;
    public final ListMultimap<Chromosome, ChrBaseRegion> PanelWithHotspots;
    public final ListMultimap<Chromosome,VariantHotspot> Hotspots;
    public final ListMultimap<Chromosome, ChrBaseRegion> HighConfidence;

    public final List<HmfTranscriptRegion> TranscriptRegions;

    public final IndexedFastaSequenceFile RefGenome;

    private final SageConfig mConfig;

    public ReferenceData(final SageConfig config)
    {
        mConfig = config;

        CoveragePanel = ArrayListMultimap.create();
        PanelWithHotspots = ArrayListMultimap.create();
        Hotspots = ArrayListMultimap.create();
        HighConfidence = ArrayListMultimap.create();

        TranscriptRegions = Lists.newArrayList();

        if(!config.AppendMode)
        {
            TranscriptRegions.addAll(config.RefGenVersion.is37() ? HmfGenePanelSupplier.allGeneList37() : HmfGenePanelSupplier.allGeneList38());
            config.Quality.populateGeneData(TranscriptRegions);
        }

        RefGenome = loadRefGenome(config.RefGenomeFile);
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

            final ListMultimap<Chromosome, ChrBaseRegion> panelWithoutHotspots = readUnnamedBedFile(mConfig.PanelBed);

            if(!mConfig.PanelBed.isEmpty())
            {
                SG_LOGGER.info("read {} panel entries from bed file: {}", panelWithoutHotspots.size(), mConfig.PanelBed);
            }

            PanelWithHotspots.putAll(panelWithHotspots(panelWithoutHotspots, Hotspots));

            if(!mConfig.HighConfidenceBed.isEmpty())
            {
                HighConfidence.putAll(readUnnamedBedFile(mConfig.HighConfidenceBed));
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

    private ListMultimap<Chromosome, ChrBaseRegion> panelWithHotspots(
            final ListMultimap<Chromosome, ChrBaseRegion> panelWithoutHotspots, final ListMultimap<Chromosome,VariantHotspot> hotspots)
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

    private static ListMultimap<Chromosome, ChrBaseRegion> readUnnamedBedFile(final String panelBed) throws IOException
    {
        final ListMultimap<Chromosome, ChrBaseRegion> panel = ArrayListMultimap.create();

        if(!panelBed.isEmpty())
        {
            SortedSetMultimap<String,GenomeRegion> bed = BEDFileLoader.fromBedFile(panelBed);

            for(String contig : bed.keySet())
            {
                if(HumanChromosome.contains(contig))
                {
                    panel.putAll(
                            HumanChromosome.fromString(contig),
                            bed.get(contig).stream().map(x -> ChrBaseRegion.from(x)).collect(Collectors.toSet()));
                }
            }
        }

        return panel;
    }
}
