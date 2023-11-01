package com.hartwig.hmftools.dnds;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.dnds.DndsCommon.DN_LOGGER;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.CodingEffect;

public class SomaticVariant
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;
    public final String Gene;
    public final boolean Biallelic;
    public final boolean Hotspot;
    public final CodingEffect WorstCodingEffect;
    public final CodingEffect CanonicalCodingEffect;
    public final int RepeatCount;

        public SomaticVariant(
            final String chromosome, final int position, final String ref, final String alt, final String gene,
            final boolean biallelic, final boolean hotspot, final CodingEffect worstCodingEffect, final CodingEffect canonicalCodingEffect,
            final int repeatCount)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Gene = gene;
        Biallelic = biallelic;
        Hotspot = hotspot;
        WorstCodingEffect = worstCodingEffect;
        CanonicalCodingEffect = canonicalCodingEffect;
        RepeatCount = repeatCount;
    }

    public static final String SAMPLE_VARIANTS_FILE = ".dnds_variants.tsv";
    public static final String COHORT_VARIANTS_FILE = "cohort_dnds_variants.tsv";

    public static String sampleDndsVariantsFilename(final String sourceDir, final String sampleId)
    {
        return sourceDir + sampleId + SAMPLE_VARIANTS_FILE;
    }

    public static String cohortDndsVariantsFilename(final String sourceDir)
    {
        return sourceDir + COHORT_VARIANTS_FILE;
    }

    public static void writeVariants(final String outputDir, final String sampleId, final List<SomaticVariant> variants)
    {
        try
        {
            String filename = sampleDndsVariantsFilename(outputDir, sampleId);
            BufferedWriter writer = createBufferedWriter(filename);

            writer.write("Chromosome\tPosition\tRef\tAlt\tGene\tBiallelic\tHotspot\tWorstCodingEffect\tCanonicalCodingEffect\tRepeatCount");
            writer.newLine();

            for(SomaticVariant variant : variants)
            {
                StringJoiner sj = new StringJoiner(TSV_DELIM);
                sj.add(variant.Chromosome).add(String.valueOf(variant.Position)).add(variant.Ref).add(variant.Alt);
                sj.add(variant.Gene).add(String.valueOf(variant.Biallelic)).add(String.valueOf(variant.Hotspot));
                sj.add(String.valueOf(variant.WorstCodingEffect)).add(String.valueOf(variant.CanonicalCodingEffect));
                sj.add(String.valueOf(variant.RepeatCount));
                writer.write(sj.toString());
                writer.newLine();
            }

            writer.close();
        }
        catch (IOException e)
        {
            DN_LOGGER.error("failed to write sample variants: {}", e.toString());
        }
    }

    // UNUSED?
    public static List<SomaticVariant> readVariants(final String filename)
    {
        List<SomaticVariant> variants = Lists.newArrayList();

        try
        {
            final List<String> lines = Files.readAllLines(new File(filename).toPath());

            final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), TSV_DELIM);
            lines.remove(0);

            int chrIndex = fieldsIndexMap.get("Chromosome");
            int posIndex = fieldsIndexMap.get("Position");
            int refIndex = fieldsIndexMap.get("Ref");
            int altIndex = fieldsIndexMap.get("Alt");
            int geneIndex = fieldsIndexMap.get("Gene");
            int biIndex = fieldsIndexMap.get("Biallelic");
            int hotspotIndex = fieldsIndexMap.get("Hotspot");
            int wceIndex = fieldsIndexMap.get("WorstCodingEffect");
            int ceIndex = fieldsIndexMap.get("CanonicalCodingEffect");
            int rcIndex = fieldsIndexMap.get("RepeatCount");
            Integer dndsImpactIndex = fieldsIndexMap.get("DndsImpact");

            for(String line : lines)
            {
                final String[] values = line.split(TSV_DELIM, -1);

                SomaticVariant variant = new SomaticVariant(
                        values[chrIndex], Integer.parseInt(values[posIndex]), values[refIndex], values[altIndex], values[geneIndex],
                        Boolean.parseBoolean(values[biIndex]), Boolean.parseBoolean(values[hotspotIndex]), CodingEffect.valueOf(values[wceIndex]),
                        CodingEffect.valueOf(values[ceIndex]), Integer.parseInt(values[rcIndex]));

                variants.add(variant);
            }
        }
        catch (IOException e)
        {
            DN_LOGGER.error("failed to read sample variants: {}", e.toString());
            return null;
        }

        return variants;
    }
}
