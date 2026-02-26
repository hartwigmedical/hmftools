package com.hartwig.hmftools.common.purple;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.segmentation.Arm;

public class ChrArmCopyNumbersFile
{
    public static final String EXTENSION = ".purple.chromosome_arm.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + EXTENSION;
    }

    public static void write(final String filename, final List<ChrArmCopyNumber> segments) throws IOException
    {
        Collection<String> lines = Lists.newArrayList();
        lines.add(tsvFileHeader());
        segments.stream().map(x -> toTSV(x)).forEach(lines::add);
        Files.write(new File(filename).toPath(), lines);
    }

    public static List<ChrArmCopyNumber> read(final String filename) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(filename).toPath());
        List<ChrArmCopyNumber> result = Lists.newArrayList();

        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), TSV_DELIM);
        lines.remove(0);

        int chrIndex = fieldsIndexMap.get(Columns.chromosome.toString());
        int armIndex = fieldsIndexMap.get(Columns.arm.toString());
        int meanCnIndex = fieldsIndexMap.get(Columns.meanCopyNumber.toString());
        int medianCnIndex = fieldsIndexMap.get(Columns.medianCopyNumber.toString());
        int minCnIndex = fieldsIndexMap.get(Columns.minCopyNumber.toString());
        int maxCnIndex = fieldsIndexMap.get(Columns.maxCopyNumber.toString());

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            ChrArmCopyNumber chrArmCopyNumber = new ChrArmCopyNumber(
                    HumanChromosome.fromString(values[chrIndex]),
                    Arm.valueOf(values[armIndex]),
                    Double.parseDouble(values[meanCnIndex]),
                    Double.parseDouble(values[medianCnIndex]),
                    Double.parseDouble(values[minCnIndex]),
                    Double.parseDouble(values[maxCnIndex]));

            result.add(chrArmCopyNumber);
        }
        return result;
    }

    private enum Columns
    {
        chromosome,
        arm,
        meanCopyNumber,
        medianCopyNumber,
        minCopyNumber,
        maxCopyNumber
    }

    private static String tsvFileHeader()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        Arrays.stream(Columns.values()).forEach(x -> sj.add(x.toString()));
        return sj.toString();
    }

    static final DecimalFormat FORMAT = new DecimalFormat("0.0000", new DecimalFormatSymbols(Locale.ENGLISH));

    private static String toTSV(final ChrArmCopyNumber chrArmCopyNumber)
    {
        return new StringJoiner(TSV_DELIM)
                .add(chrArmCopyNumber.chromosome().toString())
                .add(chrArmCopyNumber.arm().toString())
                .add(FORMAT.format(chrArmCopyNumber.meanCopyNumber()))
                .add(FORMAT.format(chrArmCopyNumber.medianCopyNumber()))
                .add(FORMAT.format(chrArmCopyNumber.minCopyNumber()))
                .add(FORMAT.format(chrArmCopyNumber.maxCopyNumber()))
                .toString();
    }

}
