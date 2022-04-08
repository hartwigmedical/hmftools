package com.hartwig.hmftools.common.cobalt;

import static com.hartwig.hmftools.common.cobalt.CobaltCommon.DELIMITER;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createGzipBufferedReader;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createGzipBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.Writer;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.Gender;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class CobaltRatioFile
{
    private static final DecimalFormat FORMAT = new DecimalFormat("#.####");

    private static final String EXTENSION = ".cobalt.ratio.tsv.gz";

    @Deprecated
    public static final String TUMOR_ONLY_REFERENCE_SAMPLE = "DIPLOID";

    @NotNull
    public static String generateFilenameForWriting(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + EXTENSION;
    }

    @NotNull
    public static String generateFilenameForReading(final String basePath, final String sample)
    {
        String filename = basePath + File.separator + sample + EXTENSION;
        if(new File(filename).exists())
        {
            return filename;
        }

        // or maybe it is not gzipped (old version)
        filename = filename.substring(0, filename.length() - 3);
        return filename;
    }

    @NotNull
    public static ListMultimap<Chromosome,CobaltRatio> read(final String filename) throws IOException
    {
        Map<Chromosome,List<CobaltRatio>> chrRatiosMap = read(filename, null, true);

        final ListMultimap<Chromosome,CobaltRatio> result = ArrayListMultimap.create();

        for(Map.Entry<Chromosome,List<CobaltRatio>> entry : chrRatiosMap.entrySet())
        {
            HumanChromosome chromosome = HumanChromosome.fromString(entry.getKey().toString());
            entry.getValue().forEach(x -> result.put(chromosome, x));
        }

        return result;
    }

    @NotNull
    public static Map<Chromosome,List<CobaltRatio>> readWithGender(final String filename, final Gender gender, boolean hasTumor)
            throws IOException
    {
        return read(filename, gender, hasTumor);
    }

    private static final String CHROMOSOME = "chromosome";
    private static final String POSITION = "position";
    private static final String REF_READ_COUNT = "referenceReadCount";
    private static final String TUMOR_READ_COUNT = "tumorReadCount";
    private static final String REF_GC_RATIO = "referenceGCRatio";
    private static final String TUMOR_GC_RATIO= "tumorGCRatio";
    private static final String REF_GC_DIP_RATIO = "referenceGCDiploidRatio";

    private static Map<Chromosome,List<CobaltRatio>> read(final String filename, final Gender gender, boolean hasTumor)
            throws IOException
    {
        Map<Chromosome,List<CobaltRatio>> chrRatiosMap = Maps.newHashMap();

        try(BufferedReader reader = filename.endsWith(".gz") ? createGzipBufferedReader(filename) : createBufferedReader(filename))
        {
            String line = reader.readLine();
            Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(line, DELIMITER);

            int chrIndex = fieldsIndexMap.get(CHROMOSOME);
            int posIndex = fieldsIndexMap.get(POSITION);
            int refReadCountIndex = fieldsIndexMap.get(REF_READ_COUNT);
            int tumorReadCountIndex = fieldsIndexMap.get(TUMOR_READ_COUNT);
            int refGcRatioIndex = fieldsIndexMap.get(REF_GC_RATIO);
            int tumorGcRatioIndex = fieldsIndexMap.get(TUMOR_GC_RATIO);
            int refGcDiplodRatioIndex = fieldsIndexMap.get(REF_GC_DIP_RATIO);

            List<CobaltRatio> ratios = null;
            String currentChromosome = "";

            while((line = reader.readLine()) != null)
            {
                String[] values = line.split(DELIMITER, -1);

                String chromosome = values[chrIndex];
                double initialRefGCRatio = Double.parseDouble(values[refGcRatioIndex]);
                double initialRefGCDiploidRatio = Double.parseDouble(values[refGcDiplodRatioIndex]);

                double refGcRatio = genderAdjustedDiploidRatio(gender, chromosome, initialRefGCRatio);
                double refGcDiploadRatio = genderAdjustedDiploidRatio(gender, chromosome, initialRefGCDiploidRatio);
                double tumorGCRatio = hasTumor ? Double.parseDouble(values[tumorGcRatioIndex]) : refGcDiploadRatio;

                CobaltRatio ratio = ImmutableCobaltRatio.builder()
                        .chromosome(chromosome)
                        .position(Integer.parseInt(values[posIndex]))
                        .referenceReadCount(Integer.parseInt(values[refReadCountIndex]))
                        .tumorReadCount(Integer.parseInt(values[tumorReadCountIndex]))
                        .tumorGCRatio(tumorGCRatio)
                        .referenceGCRatio(refGcRatio)
                        .referenceGCDiploidRatio(refGcDiploadRatio)
                        .build();

                if(!currentChromosome.equals(chromosome))
                {
                    currentChromosome = chromosome;
                    ratios = Lists.newArrayList();
                    chrRatiosMap.put(HumanChromosome.fromString(chromosome), ratios);
                }

                ratios.add(ratio);
            }
        }

        return chrRatiosMap;
    }

    public static void write(final String fileName, Collection<CobaltRatio> ratios) throws IOException
    {
        List<CobaltRatio> sorted = new ArrayList<>(ratios);
        Collections.sort(sorted);
        try(Writer writer = createGzipBufferedWriter(fileName))
        {
            for(String line : toLines(sorted))
            {
                writer.write(line + '\n');
            }
        }
    }

    private static List<String> toLines(final List<CobaltRatio> ratio)
    {
        final List<String> lines = new ArrayList<>();
        lines.add(header());
        ratio.stream().map(CobaltRatioFile::toString).forEach(lines::add);
        return lines;
    }

    private static String header()
    {
        return new StringJoiner(DELIMITER)
                .add(CHROMOSOME)
                .add(POSITION)
                .add(REF_READ_COUNT)
                .add(TUMOR_READ_COUNT)
                .add(REF_GC_RATIO)
                .add(TUMOR_GC_RATIO)
                .add(REF_GC_DIP_RATIO)
                .toString();
    }

    private static String toString(final CobaltRatio position)
    {
        return new StringJoiner(DELIMITER)
                .add(position.chromosome())
                .add(String.valueOf(position.position()))
                .add(String.valueOf(position.referenceReadCount()))
                .add(String.valueOf(position.tumorReadCount()))
                .add(FORMAT.format(position.referenceGCRatio()))
                .add(FORMAT.format(position.tumorGCRatio()))
                .add(FORMAT.format(position.referenceGCDiploidRatio()))
                .toString();
    }

    private static double genderAdjustedDiploidRatio(@Nullable final Gender gender, final String contig, double initialRatio)
    {
        if(gender == null || !HumanChromosome.contains(contig))
        {
            return initialRatio;
        }

        HumanChromosome chromosome = HumanChromosome.fromString(contig);
        if(chromosome.equals(HumanChromosome._X))
        {
            return gender.equals(Gender.FEMALE) ? 1 : 0.5;
        }

        if(chromosome.equals(HumanChromosome._Y))
        {
            return gender.equals(Gender.FEMALE) ? 0 : 0.5;
        }

        return initialRatio;
    }
}
