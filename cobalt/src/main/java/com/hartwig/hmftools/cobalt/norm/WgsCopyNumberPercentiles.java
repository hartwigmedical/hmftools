package com.hartwig.hmftools.cobalt.norm;

import static java.lang.String.format;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.BaseRegion;

public class WgsCopyNumberPercentiles
{
    private final Map<String, List<RegionPercentile>> mChrRegionsMap;

    private static final double DEFAULT_PERCENTILE = 0.5;

    private class RegionPercentile extends BaseRegion
    {
        public final double Percentile;

        public RegionPercentile(int regionStart, int regionEnd, double percentile)
        {
            super(regionStart, regionEnd);
            Percentile = percentile;
        }

        public String toString() { return format("%s: %.2f", super.toString(), Percentile); }
    }

    public WgsCopyNumberPercentiles(final String filename, final RefGenomeVersion refGenVersion)
    {
        mChrRegionsMap = Maps.newHashMap();

        loadRegions(filename, refGenVersion);
    }

    public Double getRegionPercentile(final String chromosome, final int position)
    {
        List<RegionPercentile> regionPercentiles = mChrRegionsMap.get(chromosome);

        if(regionPercentiles == null)
            return DEFAULT_PERCENTILE;

        for(RegionPercentile regionPercentile : regionPercentiles)
        {
            if(regionPercentile.containsPosition(position))
                return regionPercentile.Percentile;

            if(regionPercentile.start() > position) // in expected case there are gaps in the percentile data
                return regionPercentile.Percentile;
        }

        return DEFAULT_PERCENTILE;
    }

    private void loadRegions(final String filename, final RefGenomeVersion refGenVersion)
    {
        if(filename != null)
        {
            RefGenomeCoordinates refGenomeCoordinates = RefGenomeCoordinates.refGenomeCoordinates(refGenVersion);

            try
            {
                List<String> lines = Files.readAllLines(Paths.get(filename));

                String header = lines.get(0);
                Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
                lines.remove(0);

                int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
                int posIndex = fieldsIndexMap.get(FLD_POSITION_START);
                int percentileIndex = fieldsIndexMap.get("Percentile");

                String currentChr = "";
                int chromosomeLength = 0;
                List<RegionPercentile> regionPercentiles = null;
                RegionPercentile lastRegion = null;

                for(String line : lines)
                {
                    String[] values = line.split(TSV_DELIM, -1);

                    String chromosome = values[chrIndex];
                    int regionStart = Integer.parseInt(values[posIndex]);
                    double percentile = Double.parseDouble(values[percentileIndex]);

                    if(!chromosome.equals(currentChr))
                    {
                        if(lastRegion != null)
                            lastRegion.setEnd(chromosomeLength);

                        lastRegion = null;

                        currentChr = chromosome;
                        chromosomeLength = refGenomeCoordinates.length(chromosome);
                        regionPercentiles = Lists.newArrayList();

                        mChrRegionsMap.put(chromosome, regionPercentiles);
                    }

                    RegionPercentile region = new RegionPercentile(regionStart, regionStart, percentile);

                    if(lastRegion != null)
                        lastRegion.setEnd(regionStart - 1);

                    regionPercentiles.add(region);

                    lastRegion = region;
                }

                if(lastRegion != null)
                    lastRegion.setEnd(chromosomeLength);
            }
            catch(IOException e)
            {
                CB_LOGGER.error("failed to load WGS copy number percentiles file({}): {}", filename, e.toString());
                System.exit(1);
            }
        }
    }
}
