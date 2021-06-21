package com.hartwig.hmftools.sage;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

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
import com.hartwig.hmftools.common.genome.region.BEDFileLoader;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegionsValidation;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotFile;
import com.hartwig.hmftools.sage.config.SageConfig;

import org.apache.commons.compress.utils.Lists;

public class ReferenceData
{
    public final ListMultimap<Chromosome,NamedBed> CoveragePanel;
    public final ListMultimap<Chromosome, BaseRegion> PanelWithHotspots;
    public final ListMultimap<Chromosome,VariantHotspot> Hotspots;
    public final ListMultimap<Chromosome,BaseRegion> HighConfidence;

    private final SageConfig mConfig;

    public ReferenceData(final SageConfig config)
    {
        mConfig = config;

        CoveragePanel = ArrayListMultimap.create();
        PanelWithHotspots = ArrayListMultimap.create();
        Hotspots = ArrayListMultimap.create();
        HighConfidence = ArrayListMultimap.create();
    }

    public boolean load()
    {
        try
        {
            if(!mConfig.CoverageBed.isEmpty())
            {
                CoveragePanel.putAll(readNamedBedFile(mConfig.CoverageBed));
                SG_LOGGER.info("read {} coverage entries from bed file: {}", CoveragePanel.size(), mConfig.CoverageBed);
            }

            loadHotspots();

            final ListMultimap<Chromosome, BaseRegion> panelWithoutHotspots = readUnnamedBedFile(mConfig.PanelBed);

            if(!mConfig.PanelBed.isEmpty())
            {
                SG_LOGGER.info("read {} panel entries from bed file: {}", panelWithoutHotspots.size(), mConfig.PanelBed);
            }

            PanelWithHotspots.putAll(panelWithHotspots(panelWithoutHotspots, Hotspots));

            if(!mConfig.HighConfidenceBed.isEmpty())
            {
                HighConfidence.putAll(readUnnamedBedFile(mConfig.HighConfidenceBed));
                SG_LOGGER.info("read {} high-confidence entries from bed file: {}", HighConfidence.size(), mConfig.CoverageBed);
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

    private ListMultimap<Chromosome,BaseRegion> panelWithHotspots(
            final ListMultimap<Chromosome,BaseRegion> panelWithoutHotspots, final ListMultimap<Chromosome,VariantHotspot> hotspots)
    {
        final ListMultimap<Chromosome,BaseRegion> result = ArrayListMultimap.create();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            List<BaseRegion> regions = Lists.newArrayList();

            if(panelWithoutHotspots.containsKey(chromosome))
            {
                regions.addAll(panelWithoutHotspots.get(chromosome));
            }

            if(hotspots.containsKey(chromosome))
            {
                hotspots.get(chromosome).forEach(x -> regions.add(new BaseRegion(chromosome.toString(), (int)x.position(), (int)x.end())));
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
            int entryCount = 0;

            for(NamedBed bed : NamedBedFile.readBedFile(panelBed))
            {
                if(HumanChromosome.contains(bed.chromosome()))
                {
                    ++entryCount;
                    panel.put(HumanChromosome.fromString(bed.chromosome()), bed);
                }
            }
        }

        return panel;
    }

    private static ListMultimap<Chromosome,BaseRegion> readUnnamedBedFile(final String panelBed) throws IOException
    {
        final ListMultimap<Chromosome,BaseRegion> panel = ArrayListMultimap.create();

        if(!panelBed.isEmpty())
        {
            SortedSetMultimap<String,GenomeRegion> bed = BEDFileLoader.fromBedFile(panelBed);

            for(String contig : bed.keySet())
            {
                if(HumanChromosome.contains(contig))
                {
                    panel.putAll(
                            HumanChromosome.fromString(contig),
                            bed.get(contig).stream().map(x -> BaseRegion.from(x)).collect(Collectors.toSet()));
                }
            }
        }

        return panel;
    }
}
