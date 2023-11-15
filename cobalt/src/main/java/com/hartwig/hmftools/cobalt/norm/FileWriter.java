package com.hartwig.hmftools.cobalt.norm;

import static java.lang.String.format;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class FileWriter
{
    public static void writeNormalisationFile(
            final Map<String,List<RegionData>> chrRegionData, final RefGenomeVersion refGenomeVersion, final String outputFile)
    {
        CB_LOGGER.info("writing normalisation file: {}", outputFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("chromosome\tposition\trelativeEnrichment");
            writer.newLine();

            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                String chrStr = refGenomeVersion.versionedChromosome(chromosome.toString());

                if(!chrRegionData.containsKey(chrStr))
                    continue;

                for(RegionData regionData : chrRegionData.get(chrStr))
                {
                    writer.write(format("%s\t%d\t%.4f",
                            chrStr, regionData.Position,
                            regionData.relativeEnrichment() > 0 ? regionData.relativeEnrichment() : Double.NaN));

                    writer.newLine();
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            CB_LOGGER.error("failed to write normalisation file: {}", e.toString());
            System.exit(1);
        }
    }

    public static void writeDetailedFile(
            final Map<String,List<RegionData>> chrRegionData, final List<String> sampleIds, final String outputFile)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("SampleId\tChromosome\tPosition\tGcBucket\tMappability\tGcRatioPanel\tReadDepth\tAdjGcRatio");
            writer.newLine();

            for(Map.Entry<String, List<RegionData>> entry : chrRegionData.entrySet())
            {
                String chromosome = entry.getKey();

                List<RegionData> regions = entry.getValue();

                for(RegionData regionData : regions)
                {
                    for(int i = 0; i < sampleIds.size(); ++i)
                    {
                        String sampleId = sampleIds.get(i);
                        SampleRegionData sampleRegionData = regionData.getSampleData(i);

                        writer.write(format("%s\t%s\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f",
                                sampleId, chromosome, regionData.Position, regionData.gcBucket(), regionData.mappability(),
                                sampleRegionData.GcRatioPanel, sampleRegionData.ReadDepth, sampleRegionData.adjustedGcRatio()));
                        writer.newLine();
                    }
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            CB_LOGGER.error("failed to write normalisation file: {}", e.toString());
            System.exit(1);
        }
    }

}
