package com.hartwig.hmftools.common.amber;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.*;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.io.Writer;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.primitives.Doubles;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.file.FileReaderUtils;

public final class AmberBAFFile
{
    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");

    private static final String AMBER_EXTENSION = ".amber.baf.tsv.gz";
    private static final String AMBER_EXTENSION_OLD = ".amber.baf.tsv";

    public static String generateAmberFilenameForWriting(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + AMBER_EXTENSION;
    }

    public static String generateAmberFilenameForReading(final String basePath, final String sample)
    {
        String filename = checkAddDirSeparator(basePath) + sample + AMBER_EXTENSION;

        if(Files.exists(Paths.get(filename)))
            return filename;

        return checkAddDirSeparator(basePath) + sample + AMBER_EXTENSION_OLD;
    }

    private static final String CHROMOSOME = "chromosome";
    private static final String POSITION = "position";
    private static final String TUMOR_BAF = "tumorBAF";
    private static final String TUMOR_MOD_BAF = "tumorModifiedBAF";
    private static final String TUMOR_DEPTH = "tumorDepth";
    private static final String NORM_BAF = "normalBAF";
    private static final String NORM_MOD_BAF = "normalModifiedBAF";
    private static final String NORM_DEPTH = "normalDepth";

    public static Multimap<Chromosome,AmberBAF> read(final String fileName, boolean hasTumor) throws IOException
    {
        ListMultimap<Chromosome,AmberBAF> chrBafMap = ArrayListMultimap.create();

        try(BufferedReader reader = fileName.endsWith(".gz") ? createGzipBufferedReader(fileName) : createBufferedReader(fileName))
        {
            String line = reader.readLine();
            Map<String, Integer> fieldsIndexMap = FileReaderUtils.createFieldsIndexMap(line, TSV_DELIM);

            int chrIndex = fieldsIndexMap.get(CHROMOSOME);
            int posIndex = fieldsIndexMap.get(POSITION);
            int tumorBafIndex = fieldsIndexMap.get(TUMOR_BAF);
            int tumorDepthIndex = fieldsIndexMap.get(TUMOR_DEPTH);
            int normBafIndex = fieldsIndexMap.get(NORM_BAF);
            int normDepthIndex = fieldsIndexMap.get(NORM_DEPTH);

            while((line = reader.readLine()) != null)
            {
                String[] values = line.split(TSV_DELIM, -1);

                String chromosome = values[chrIndex];

                double tumorBAF = hasTumor ? Double.parseDouble(values[tumorBafIndex]) : 0.5;

                AmberBAF amberBAF = new AmberBAF(
                        chromosome, Integer.parseInt(values[posIndex]), tumorBAF, Integer.parseInt(values[tumorDepthIndex]),
                        Double.parseDouble(values[normBafIndex]), Integer.parseInt(values[normDepthIndex]));

                chrBafMap.put(HumanChromosome.fromString(chromosome), amberBAF);
            }
        }

        return chrBafMap;
    }

    public static void write(final String filename, final List<AmberBAF> bafs) throws IOException
    {
        try (Writer writer = createGzipBufferedWriter(filename))
        {
            for (String line : toLines(bafs))
            {
                writer.write(line + '\n');
            }
        }
    }

    private static List<String> toLines(final List<AmberBAF> bafs)
    {
        final List<String> lines = new ArrayList<>();
        lines.add(header());
        bafs.stream().map(AmberBAFFile::toString).forEach(lines::add);
        return lines;
    }

    private static String header()
    {
        return new StringJoiner(TSV_DELIM, "", "")
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
        return new StringJoiner(TSV_DELIM).add(baf.Chromosome)
                .add(String.valueOf(baf.Position))
                .add(FORMAT.format(baf.tumorBAF()))
                .add(FORMAT.format(baf.tumorModifiedBAF()))
                .add(String.valueOf(baf.tumorDepth()))
                .add(Doubles.isFinite(baf.normalBAF()) ? FORMAT.format(baf.normalBAF()) : "0")
                .add(Doubles.isFinite(baf.normalModifiedBAF()) ? FORMAT.format(baf.normalModifiedBAF()) : "0")
                .add(String.valueOf(baf.normalDepth()))
                .toString();
    }
}
