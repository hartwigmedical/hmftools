package com.hartwig.hmftools.common.cobalt;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ZIP_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createGzipBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfileFactory;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class CobaltRatioFile
{
    public enum Column
    {
        chromosome,
        position,
        referenceReadDepth,
        tumorReadDepth,
        referenceGCRatio,
        tumorGCRatio,
        referenceGCDiploidRatio,
        referenceGCContent,
        tumorGCContent
    }

    public static final DecimalFormat FORMAT = new DecimalFormat("#.####", new DecimalFormatSymbols(Locale.ENGLISH));

    private static final String EXTENSION = ".cobalt.ratio.tsv.gz";
    private static final String EXTENSION_UNZIPPED = ".cobalt.ratio.tsv";

    // old column names for backwards compatibility
    private static final String COL_REF_READ_COUNT  = "referenceReadCount";
    private static final String COL_TUMOR_READ_COUNT  = "tumorReadCount";

    @Deprecated
    public static final String TUMOR_ONLY_REFERENCE_SAMPLE = "DIPLOID";

    @NotNull
    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + EXTENSION;
    }

    @NotNull
    public static String generateFilenameUnzipped(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + EXTENSION_UNZIPPED;
    }

    @NotNull
    public static String generateFilenameForReading(final String basePath, final String sample)
    {
        // some old samples have unzipped ratio files, so check for these
        String filename = generateFilenameUnzipped(basePath, sample);

        if(Files.exists(Paths.get(filename)))
            return filename;

        String unzippedFile = filename.replaceAll(ZIP_EXTENSION, "");

        if(Files.exists(Paths.get(unzippedFile)))
            return unzippedFile;

        return filename;
    }

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

    public static Map<Chromosome,List<CobaltRatio>> readWithGender(final String filename, final Gender gender, boolean hasTumor)
            throws IOException
    {
        return read(filename, gender, hasTumor);
    }

    private static final int DEFAULT_READ_LENGTH = 151;
    public static final double READ_DEPTH_INVALID = -1;

    private static double convertReadCount(final double readCount)
    {
        if(readCount <= 0)
            return readCount;

        return readCount * DEFAULT_READ_LENGTH /GCProfileFactory.WINDOW_SIZE;
    }

    private static Map<Chromosome,List<CobaltRatio>> read(final String filename, final Gender gender, boolean hasTumor)
    {
        Map<Chromosome,List<CobaltRatio>> chrRatiosMap = new HashMap<>();

        try(DelimFileReader reader = new DelimFileReader(filename))
        {
            List<CobaltRatio> ratios = null;
            String currentChromosome = null;

            int chrIndex = reader.getColumnIndex(Column.chromosome);
            int posIndex = reader.getColumnIndex(Column.position);

            int refGcRatioIndex = reader.getColumnIndex(Column.referenceGCRatio);
            int tumorGcRatioIndex = reader.getColumnIndex(Column.tumorGCRatio);
            int refGcDiplodRatioIndex = reader.getColumnIndex(Column.referenceGCDiploidRatio);

            // 1.15 field updates
            // changes: referenceReadCount -> referenceReadDepth, tumorReadCount -> tumorReadDepth
            // added: referenceGCContent, tumorGCContent

            // v1.16 onwards
            Integer refReadDepthIndex = reader.getColumnIndex(Column.referenceReadDepth);
            Integer tumorReadDepthIndex = reader.getColumnIndex(Column.tumorReadDepth);
            Integer refGcContentIndex = reader.getColumnIndex(Column.referenceGCContent);
            Integer tumorGcContentIndex = reader.getColumnIndex(Column.tumorGCContent);

            // v1.15 backwards compatibility with conversion below
            Integer refReadCountIndex = reader.getColumnIndex(COL_REF_READ_COUNT);
            Integer tumorReadCountIndex = reader.getColumnIndex(COL_TUMOR_READ_COUNT);

            boolean useReadCount = refReadCountIndex != null && tumorReadCountIndex != null;
            boolean hasGcContent = refGcContentIndex != null && tumorGcContentIndex != null;

            for(DelimFileReader.Row row : reader)
            {
                String chromosome = row.get(chrIndex);

                if(currentChromosome == null || !currentChromosome.equals(chromosome))
                {
                    currentChromosome = chromosome;
                    ratios = new ArrayList<>();
                    chrRatiosMap.put(HumanChromosome.fromString(chromosome), ratios);
                }
                else
                {
                    // use the same String object for better performance
                    chromosome = currentChromosome;
                }

                double refReadDepth = useReadCount ? convertReadCount(row.getDouble(refReadCountIndex)) : row.getDouble(refReadDepthIndex);

                double initialRefGCRatio = row.getDouble(refGcRatioIndex);
                double initialRefGCDiploidRatio = row.getDouble(refGcDiplodRatioIndex);

                if(refReadDepth == READ_DEPTH_INVALID)
                {
                    // revert to a default ref ratio where no information is available (ie in tumor/panel only)
                    initialRefGCRatio = 1;
                    initialRefGCDiploidRatio = 1;
                }

                double refGcRatio = genderAdjustedDiploidRatio(gender, chromosome, initialRefGCRatio);
                double refGcDiploadRatio = genderAdjustedDiploidRatio(gender, chromosome, initialRefGCDiploidRatio);
                double tumorGCRatio = hasTumor ? row.getDouble(tumorGcRatioIndex) : refGcDiploadRatio;

                double tumorReadDepth = useReadCount ? convertReadCount(row.getDouble(tumorReadCountIndex)) : row.getDouble(tumorReadDepthIndex);

                double refGcContent = hasGcContent ? row.getDouble(refGcContentIndex) : 0;
                double tumorGcContent = hasGcContent ? row.getDouble(tumorGcContentIndex) : 0;

                CobaltRatio ratio = new CobaltRatio(
                        chromosome,
                        row.getInt(posIndex),
                        refReadDepth,
                        refGcRatio,
                        refGcContent,
                        refGcDiploadRatio,
                        tumorReadDepth,
                        tumorGCRatio,
                        tumorGcContent);

                ratios.add(ratio);
            }
        }

        return chrRatiosMap;
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

    public static void write(final String fileName, Collection<CobaltRatio> ratios) throws IOException
    {
        List<CobaltRatio> sorted = new ArrayList<>(ratios);
        Collections.sort(sorted);

        try(BufferedWriter writer = createGzipBufferedWriter(fileName))
        {
            DelimFileWriter.write(writer, Column.values(), sorted,
                (ratio, row) -> {
                    row.set(Column.chromosome, ratio.chromosome());
                    row.set(Column.position, ratio.position());
                    row.set(Column.referenceReadDepth, ratio.referenceReadDepth(), FORMAT);
                    row.set(Column.tumorReadDepth, ratio.tumorReadDepth(), FORMAT);
                    row.set(Column.referenceGCRatio, ratio.referenceGCRatio(), FORMAT);
                    row.set(Column.tumorGCRatio, ratio.tumorGCRatio(), FORMAT);
                    row.set(Column.referenceGCDiploidRatio, ratio.referenceGCDiploidRatio(), FORMAT);
                    row.set(Column.referenceGCContent, ratio.referenceGcContent(), FORMAT);
                    row.set(Column.tumorGCContent, ratio.tumorGcContent(), FORMAT);
                });
        }
    }
}
