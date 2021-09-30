package com.hartwig.hmftools.common.purple.cnchromosome;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.segment.ChromosomeArm;

import org.jetbrains.annotations.NotNull;

public final class CnPerChromosomeArmFile
{

    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");
    private static final String COPYNUMBER_PER_CHROMOSOME_ARM_EXTENSION = ".cnv.chromosomearm.somatic.tsv";
    private static final String DELIMITER = "\t";

    public static String generateFilename(@NotNull String basePath, @NotNull String sample)
    {
        return basePath + File.separator + sample + COPYNUMBER_PER_CHROMOSOME_ARM_EXTENSION;
    }

    public static void write(@NotNull String filename, @NotNull List<CnPerChromosomeArmData> cnPerChromosomeArm) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(cnPerChromosomeArm));
    }

    private static List<String> toLines(@NotNull List<CnPerChromosomeArmData> cnPerChromosomeArm)
    {
        List<String> lines = Lists.newArrayList();
        lines.add(header());

        for(CnPerChromosomeArmData chromosomeArm : cnPerChromosomeArm)
        {
            lines.add(new StringJoiner(DELIMITER).add(chromosomeArm.chromosome().toString())
                    .add(chromosomeArm.chromosomeArm().name())
                    .add(FORMAT.format(chromosomeArm.copyNumber()))
                    .toString());
        }
        return lines;
    }

    private static String header()
    {
        return new StringJoiner(DELIMITER)
                .add("chromosome")
                .add("chromosomeArm")
                .add("copyNumber")
                .toString();
    }

    public static List<CnPerChromosomeArmData> fromLines(final String filename) throws IOException
    {
        List<CnPerChromosomeArmData> armDataList = Lists.newArrayList();
        final List<String> lines = Files.readAllLines(Paths.get(filename));

        String header = lines.get(0);
        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIMITER);
        lines.remove(0);

        for(String line : lines)
        {
            String[] values = line.split(DELIMITER, -1);

            armDataList.add(ImmutableCnPerChromosomeArmData.builder()
                    .chromosome(HumanChromosome.fromString(values[fieldsIndexMap.get("chromosome")]))
                    .chromosomeArm(ChromosomeArm.valueOf(values[fieldsIndexMap.get("chromosomeArm")]))
                    .copyNumber(Double.parseDouble(values[fieldsIndexMap.get("copyNumber")]))
                    .build());
        }

        return armDataList;
    }

}
