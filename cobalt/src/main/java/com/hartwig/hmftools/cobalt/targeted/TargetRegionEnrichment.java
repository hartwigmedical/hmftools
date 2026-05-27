package com.hartwig.hmftools.cobalt.targeted;

import static java.lang.String.format;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.StringJoiner;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class TargetRegionEnrichment
{
    public final Chromosome Chromosome;
    public final int Position;
    public final double Enrichment;

    public TargetRegionEnrichment(final Chromosome chromosome, final int position, final double enrichment)
    {
        Chromosome = chromosome;
        Position = position;
        Enrichment = enrichment;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final TargetRegionEnrichment that = (TargetRegionEnrichment) o;
        return Position == that.Position && Objects.equals(Chromosome, that.Chromosome);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Chromosome, Position);
    }

    @Override
    public String toString()
    {
        return format("%s:%d enrichment(%.4f)", Chromosome, Position, Enrichment);
    }

    private enum Columns
    {
        Chromosome,
        Position,
        RelativeEnrichment;
    }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        Arrays.stream(Columns.values()).forEach(x -> sj.add(x.toString()));
        return sj.toString();
    }

    public static String toTsv(final String chromosome, final int position, final double enrichment)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(chromosome);
        sj.add(String.valueOf(position));
        sj.add(String.format("%.4f", enrichment));
        return sj.toString();
    }

    public static ListMultimap<Chromosome,TargetRegionEnrichment> loadEnrichmentFile(final String filename)
    {
        ListMultimap<Chromosome, TargetRegionEnrichment> result = ArrayListMultimap.create();

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));

            String header = lines.get(0);
            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
            lines.remove(0);

            int chrIndex = fieldsIndexMap.getOrDefault(Columns.Chromosome.toString(), 0);
            int posIndex = fieldsIndexMap.getOrDefault(Columns.Position.toString(), 1);
            int enrichmentIndex = fieldsIndexMap.getOrDefault(Columns.RelativeEnrichment.toString(), 2);

            for(String line : lines)
            {
                String[] values = line.split(TSV_DELIM, -1);

                HumanChromosome chromosome = HumanChromosome.fromString(values[chrIndex]);

                result.put(chromosome, new TargetRegionEnrichment(
                        chromosome, Integer.parseInt(values[posIndex]), Double.parseDouble(values[enrichmentIndex])));
            }

            CB_LOGGER.info("loaded {} target-panel norm regions from file({})", result.size(), filename);
            return result;
        }
        catch(IOException e)
        {
            CB_LOGGER.error("failed to load target-panel norm regions file({}): {}", filename, e.toString());
            return null;
        }
    }

}
