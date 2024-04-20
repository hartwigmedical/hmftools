package com.hartwig.hmftools.cobalt.norm;

import static java.lang.Math.floor;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.norm.NormConstants.REGION_SIZE;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class DataLoader
{
    public static void addTargetRegions(final List<ChrBaseRegion> bedRegions, final Map<String,List<RegionData>> chrRegionData)
    {
        String currentChromosome = "";
        List<RegionData> regions = null;

        for(ChrBaseRegion region : bedRegions)
        {
            if(!region.chromosome().equals(currentChromosome))
            {
                currentChromosome = region.chromosome();
                regions = Lists.newArrayList();
                chrRegionData.put(region.chromosome(), regions);
            }

            int startPosition = (int)(floor(region.start()/(double)REGION_SIZE) * REGION_SIZE + 1);
            int endPosition = startPosition + REGION_SIZE - 1;

            addRegion(regions, startPosition);

            while(endPosition < region.end())
            {
                startPosition += REGION_SIZE;
                addRegion(regions, startPosition);
                endPosition = startPosition + REGION_SIZE - 1;
            }
        }
    }

    private static void addRegion(final List<RegionData> regions, int position)
    {
        RegionData prevRegion = !regions.isEmpty() ? regions.get(regions.size() - 1) : null;
        if(prevRegion != null && prevRegion.Position == position)
            return;

        regions.add(new RegionData(position));
    }

    public static void addCobaltSampleData(
            final Gender amberGender, final String cobaltPanelFilename, final String cobaltWgsFilename,
            final Map<String,List<RegionData>> chrRegionData)
    {
        try
        {
            CB_LOGGER.info("reading Cobalt ratios from {}", cobaltPanelFilename);

            Map<Chromosome,List<CobaltRatio>> chrPanelRatios = CobaltRatioFile.readWithGender(
                    cobaltPanelFilename, amberGender, true);

            Map<Chromosome,List<CobaltRatio>> chrWgsRatios = null;

            if(!cobaltWgsFilename.isEmpty())
            {
                CB_LOGGER.info("reading Cobalt WGS ratios from {}", cobaltWgsFilename);

                chrWgsRatios = CobaltRatioFile.readWithGender(cobaltWgsFilename, amberGender, true);
            }

            for(Map.Entry<String,List<RegionData>> entry : chrRegionData.entrySet())
            {
                String chrStr = entry.getKey();
                List<RegionData> regions = entry.getValue();

                HumanChromosome chromosome = HumanChromosome.fromString(chrStr);
                List<CobaltRatio> cobaltRatios = chrPanelRatios.get(chromosome);

                if(cobaltRatios == null)
                    continue;

                double defaultWgsGcRatio = wgsGcRatio(amberGender, chromosome);
                List<CobaltRatio> cobaltWgsRatios = chrWgsRatios != null ? chrWgsRatios.get(chromosome) : Collections.emptyList();

                int cobaltIndex = 0;
                int cobaltWgsIndex = 0;

                for(RegionData region : regions)
                {
                    CobaltRatio cobaltRatio = null;

                    while(true)
                    {
                        if(cobaltIndex >= cobaltRatios.size())
                            break;

                        cobaltRatio = cobaltRatios.get(cobaltIndex);

                        if(cobaltRatio.position() == region.Position)
                            break;
                        else if(cobaltRatio.position() > region.Position)
                            break;

                        ++cobaltIndex;
                    }

                    // likewise find the matching WGS ratio if available
                    CobaltRatio cobaltWgsRatio = null;

                    while(true)
                    {
                        if(cobaltWgsIndex >= cobaltWgsRatios.size())
                            break;

                        cobaltWgsRatio = cobaltWgsRatios.get(cobaltWgsIndex);

                        if(cobaltWgsRatio.position() == region.Position)
                            break;
                        else if(cobaltWgsRatio.position() > region.Position)
                            break;

                        ++cobaltWgsIndex;
                    }

                    double wgsGcRatio = cobaltWgsRatio != null && cobaltWgsRatio.position() == region.Position ?
                            cobaltWgsRatio.tumorGCRatio() : defaultWgsGcRatio;

                    if(cobaltRatio != null && cobaltRatio.position() == region.Position)
                    {
                        region.addSampleRegionData(new SampleRegionData(cobaltRatio.tumorReadDepth(), cobaltRatio.tumorGCRatio(), wgsGcRatio));
                    }
                    else
                    {
                        region.addSampleRegionData(new SampleRegionData(0, 0, wgsGcRatio));
                    }
                }
            }
        }
        catch(IOException e)
        {
            CB_LOGGER.error("sample({}) failed to read Cobalt data: {}", cobaltPanelFilename, e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }

    private static double wgsGcRatio(final Gender amberGender, final HumanChromosome chromosome)
    {
        if(chromosome == HumanChromosome._X)
            return amberGender == Gender.FEMALE ? 1 : 0.5;
        else if(chromosome == HumanChromosome._Y)
            return amberGender == Gender.FEMALE ? 0 : 0.5;
        else
            return 1;
    }
}
