package com.hartwig.hmftools.common.amber;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.zip.GZIPOutputStream;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.primitives.Doubles;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;

public final class AmberBAFFile
{
    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");

    private static final String DELIMITER = "\t";
    private static final String AMBER_EXTENSION = ".amber.baf.tsv.gz";
    private static final String AMBER_EXTENSION_OLD = ".amber.baf.tsv";

    @NotNull
    public static String generateAmberFilenameForWriting(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + AMBER_EXTENSION;
    }

    @NotNull
    public static String generateAmberFilenameForReading(final String basePath, final String sample)
    {
        String filename = basePath + File.separator + sample + AMBER_EXTENSION;

        if(Files.exists(Paths.get(filename)))
            return filename;

        return basePath + File.separator + sample + AMBER_EXTENSION_OLD;
    }

    private static final String CHROMOSOME = "chromosome";
    private static final String POSITION = "position";
    private static final String TUMOR_BAF = "tumorBAF";
    private static final String TUMOR_MOD_BAF = "tumorModifiedBAF";
    private static final String TUMOR_DEPTH = "tumorDepth";
    private static final String NORM_BAF = "normalBAF";
    private static final String NORM_MOD_BAF = "normalModifiedBAF";
    private static final String NORM_DEPTH = "normalDepth";

    @NotNull
    public static Multimap<Chromosome,AmberBAF> read(final String fileName, boolean hasTumor) throws IOException
    {
        ListMultimap<Chromosome,AmberBAF> chrBafMap = ArrayListMultimap.create();

        try(BufferedReader reader = createBufferedReader(fileName))
        {

            String line = reader.readLine();
            Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(line, DELIMITER);

            int chrIndex = fieldsIndexMap.get(CHROMOSOME);
            int posIndex = fieldsIndexMap.get(POSITION);
            int tumorBafIndex = fieldsIndexMap.get(TUMOR_BAF);
            int tumorDepthIndex = fieldsIndexMap.get(TUMOR_DEPTH);
            int normBafIndex = fieldsIndexMap.get(NORM_BAF);
            int normDepthIndex = fieldsIndexMap.get(NORM_DEPTH);

            while((line = reader.readLine()) != null)
            {
                String[] values = line.split(DELIMITER, -1);

                String chromosome = values[chrIndex];

                double tumorBAF = hasTumor ? Double.parseDouble(values[tumorBafIndex]) : 0.5;

                AmberBAF amberBAF = ImmutableAmberBAF.builder()
                        .chromosome(chromosome)
                        .position(Integer.parseInt(values[posIndex]))
                        .tumorBAF(tumorBAF)
                        .tumorDepth(Integer.parseInt(values[tumorDepthIndex]))
                        .normalBAF(Double.parseDouble(values[normBafIndex]))
                        .normalDepth(Integer.parseInt(values[normDepthIndex]))
                        .build();

                chrBafMap.put(HumanChromosome.fromString(chromosome), amberBAF);
            }
        }

        return chrBafMap;
    }

    public static void write(final String filename, final List<AmberBAF> bafs) throws IOException
    {
        OutputStream outputStream = new FileOutputStream(filename);
        if (filename.endsWith(".gz"))
            outputStream = new GZIPOutputStream(outputStream);
        try (Writer writer = new OutputStreamWriter(outputStream, StandardCharsets.UTF_8))
        {
            for (String line : toLines(bafs))
            {
                writer.write(line + '\n');
            }
        }
    }

    @NotNull
    private static List<String> toLines(final List<AmberBAF> bafs)
    {
        final List<String> lines = new ArrayList<>();
        lines.add(header());
        bafs.stream().map(AmberBAFFile::toString).forEach(lines::add);
        return lines;
    }

    private static String header()
    {
        return new StringJoiner(DELIMITER, "", "")
                .add(CHROMOSOME)
                .add(POSITION)
                .add(TUMOR_BAF)
                .add(TUMOR_MOD_BAF)
                .add(TUMOR_DEPTH)
                .add(NORM_BAF)
                .add(NORM_MOD_BAF)
                .add(NORM_DEPTH)
                .toString();
    }

    private static String toString(final AmberBAF baf)
    {
        return new StringJoiner(DELIMITER).add(baf.chromosome())
                .add(String.valueOf(baf.position()))
                .add(FORMAT.format(baf.tumorBAF()))
                .add(FORMAT.format(baf.tumorModifiedBAF()))
                .add(String.valueOf(baf.tumorDepth()))
                .add(Doubles.isFinite(baf.normalBAF()) ? FORMAT.format(baf.normalBAF()) : "0")
                .add(Doubles.isFinite(baf.normalModifiedBAF()) ? FORMAT.format(baf.normalModifiedBAF()) : "0")
                .add(String.valueOf(baf.normalDepth()))
                .toString();
    }
}
