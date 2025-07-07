package com.hartwig.hmftools.cobalt.norm;

import static java.lang.String.format;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

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

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add("SampleId").add("Chromosome").add("Position").add("ProfileGcBucket").add("Mappability");
            sj.add("PanelGcContent").add("ReadDepth").add("PanelGcRatio").add("AdjGcRatio");
            writer.write(sj.toString());
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

                        sj = new StringJoiner(TSV_DELIM);
                        sj.add(sampleId);
                        sj.add(chromosome);
                        sj.add(String.valueOf(regionData.Position));
                        sj.add(String.valueOf(regionData.profileGcBucket()));
                        sj.add(format("%.3f", regionData.mappability()));
                        sj.add(format("%.3f", sampleRegionData.PanelGcContent));
                        sj.add(format("%.3f", sampleRegionData.ReadDepth));
                        sj.add(format("%.3f", sampleRegionData.PanelGcRatio));
                        sj.add(format("%.3f", sampleRegionData.adjustedGcRatio()));

                        writer.write(sj.toString());
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
