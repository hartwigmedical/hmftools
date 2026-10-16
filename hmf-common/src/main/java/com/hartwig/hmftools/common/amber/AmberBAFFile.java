package com.hartwig.hmftools.common.amber;

import static java.lang.String.valueOf;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.*;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.io.Writer;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
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

    public static String generateAmberFilenameForWriting(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + AMBER_EXTENSION;
    }

    public static String generateAmberFilenameForReading(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + AMBER_EXTENSION;
    }

    private enum Column
    {
        chromosome,
        position,
        tumorBAF,
        tumorModifiedBAF,
        tumorDepth,
        normalBAF,
        normalModifiedBAF,
        normalDepth,
        lohRecovered;
    }

    public static ListMultimap<Chromosome, AmberBAF> read(final String fileName, boolean hasTumor) throws IOException
    {
        ListMultimap<Chromosome, AmberBAF> chrBafMap = ArrayListMultimap.create();

        try(BufferedReader reader = fileName.endsWith(".gz") ? createGzipBufferedReader(fileName) : createBufferedReader(fileName))
        {
            String line = reader.readLine();
            Map<String, Integer> fieldsIndexMap = FileReaderUtils.createFieldsIndexMap(line, TSV_DELIM);

            int chrIndex = fieldsIndexMap.get(Column.chromosome.toString());
            int posIndex = fieldsIndexMap.get(Column.position.toString());
            int tumorBafIndex = fieldsIndexMap.get(Column.tumorBAF.toString());
            int tumorDepthIndex = fieldsIndexMap.get(Column.tumorDepth.toString());
            int normBafIndex = fieldsIndexMap.get(Column.normalBAF.toString());
            int normDepthIndex = fieldsIndexMap.get(Column.normalDepth.toString());

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
        try(Writer writer = createGzipBufferedWriter(filename))
        {
            for(String line : toLines(bafs))
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
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        Arrays.stream(Column.values()).forEach(x -> sj.add(x.toString()));
        return sj.toString();
    }

    private static String toString(final AmberBAF baf)
    {
        return new StringJoiner(TSV_DELIM)
                .add(baf.Chromosome)
                .add(valueOf(baf.Position))
                .add(FORMAT.format(baf.tumorBAF()))
                .add(FORMAT.format(baf.tumorModifiedBAF()))
                .add(valueOf(baf.tumorDepth()))
                .add(Doubles.isFinite(baf.normalBAF()) ? FORMAT.format(baf.normalBAF()) : "0")
                .add(Doubles.isFinite(baf.normalModifiedBAF()) ? FORMAT.format(baf.normalModifiedBAF()) : "0")
                .add(valueOf(baf.normalDepth()))
                .add(valueOf(baf.LohRecovered))
                .toString();
    }
}
