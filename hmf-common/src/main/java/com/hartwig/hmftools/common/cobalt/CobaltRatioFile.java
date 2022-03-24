package com.hartwig.hmftools.common.cobalt;

import static com.hartwig.hmftools.common.cobalt.CobaltCommon.DELIMITER;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.zip.GZIPOutputStream;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class CobaltRatioFile
{
    private static final DecimalFormat FORMAT = new DecimalFormat("#.####");

    private static final String EXTENSION = ".cobalt.ratio.tsv";

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

        // or maybe it is gzipped
        filename += ".gz";
        return filename;
    }

    @NotNull
    public static ListMultimap<Chromosome, CobaltRatio> read(final String filename) throws IOException
    {
        return read(filename, null);
    }

    @NotNull
    public static ListMultimap<Chromosome, CobaltRatio> readTumorOnly(final String filename, final Gender gender)
            throws IOException
    {
        return read(filename, gender);
    }

    private static final String CHROMOSOME = "chromosome";
    private static final String POSITION = "position";
    private static final String REF_READ_COUNT = "referenceReadCount";
    private static final String TUMOR_READ_COUNT = "tumorReadCount";
    private static final String REF_GC_RATIO = "referenceGCRatio";
    private static final String TUMOR_GC_RATIO= "tumorGCRatio";
    private static final String REF_GC_DIP_RATIO = "referenceGCDiploidRatio";

    private static ListMultimap<Chromosome,CobaltRatio> read(final String filename, final Gender gender)
            throws IOException
    {
        ListMultimap<Chromosome,CobaltRatio> chrRatiosMap = ArrayListMultimap.create();

        BufferedReader reader = createBufferedReader(filename);

        String line = reader.readLine();
        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(line, DELIMITER);

        int chrIndex = fieldsIndexMap.get(CHROMOSOME);
        int posIndex = fieldsIndexMap.get(POSITION);
        int refReadCountIndex = fieldsIndexMap.get(REF_READ_COUNT);
        int tumorReadCountIndex = fieldsIndexMap.get(TUMOR_READ_COUNT);
        int refGcRatioIndex = fieldsIndexMap.get(REF_GC_RATIO);
        int tumorGcRatioIndex = fieldsIndexMap.get(TUMOR_GC_RATIO);
        int refGcDiplodRatioIndex = fieldsIndexMap.get(REF_GC_DIP_RATIO);

        while((line = reader.readLine()) != null)
        {
            String[] values = line.split(DELIMITER, -1);

            String chromosome = values[chrIndex];
            double initialRefGCRatio = Double.parseDouble(values[refGcRatioIndex]);
            double initialRefGCDiploidRatio = Double.parseDouble(values[refGcDiplodRatioIndex]);

            CobaltRatio ratio = ImmutableCobaltRatio.builder()
                    .chromosome(chromosome)
                    .position(Integer.parseInt(values[posIndex]))
                    .referenceReadCount(Integer.parseInt(values[refReadCountIndex]))
                    .tumorReadCount(Integer.parseInt(values[tumorReadCountIndex].trim()))
                    .referenceGCRatio(genderAdjustedDiploidRatio(gender, chromosome, initialRefGCRatio))
                    .tumorGCRatio(Double.parseDouble(values[tumorGcRatioIndex].trim()))
                    .referenceGCDiploidRatio(genderAdjustedDiploidRatio(gender, chromosome, initialRefGCDiploidRatio))
                    .build();

            chrRatiosMap.put(HumanChromosome.fromString(chromosome), ratio);
        }

        return chrRatiosMap;
    }

    public static void write(final String fileName, Multimap<Chromosome, CobaltRatio> ratios) throws IOException
    {
        List<CobaltRatio> sorted = new ArrayList<>(ratios.values());
        Collections.sort(sorted);
        write(fileName, sorted);
    }

    private static void write(final String fileName, List<CobaltRatio> ratios) throws IOException
    {
        OutputStream outputStream = new FileOutputStream(fileName);
        if(fileName.endsWith(".gz"))
        {
            outputStream = new GZIPOutputStream(outputStream);
        }
        try(Writer writer = new OutputStreamWriter(outputStream, StandardCharsets.UTF_8))
        {
            for(String line : toLines(ratios))
            {
                writer.write(line + '\n');
            }
        }
    }

    private static List<String> toLines(final List<CobaltRatio> ratio)
    {
        final List<String> lines = Lists.newArrayList();
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
                .add(String.valueOf(FORMAT.format(position.referenceGCRatio())))
                .add(String.valueOf(FORMAT.format(position.tumorGCRatio())))
                .add(String.valueOf(FORMAT.format(position.referenceGCDiploidRatio())))
                .toString();
    }

    private static double genderAdjustedDiploidRatio(@Nullable final Gender gender, final String contig, double initialRatio)
    {
        if(gender == null || Doubles.lessOrEqual(initialRatio, 0) || !HumanChromosome.contains(contig))
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
