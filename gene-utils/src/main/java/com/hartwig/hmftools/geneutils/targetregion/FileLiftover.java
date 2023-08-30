package com.hartwig.hmftools.geneutils.targetregion;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache.UNMAPPED_POSITION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.file.FileDelimiters;
import com.hartwig.hmftools.common.utils.file.FileReaderUtils;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class FileLiftover
{
    private static final String INPUT_FILE = "input_file";
    private static final String OUTPUT_FILE = "output_file";
    private static final String SOURCE_REF_GENOME_VERSION = "source_ref_genome_version";

    public FileLiftover(final ConfigBuilder configBuilder)
    {
        String inputFile = configBuilder.getValue(INPUT_FILE);
        String outputFile = configBuilder.getValue(OUTPUT_FILE);

        RefGenomeVersion sourceVersion = RefGenomeVersion.from(configBuilder.getValue(SOURCE_REF_GENOME_VERSION));
        RefGenomeVersion destVersion = sourceVersion.is37() ? V38 : V37;

        GenomeLiftoverCache liftoverCache = new GenomeLiftoverCache(true, destVersion == V38);

        GU_LOGGER.info("lifting over file({}) to ref genome version({})", inputFile, destVersion);

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(inputFile));
            String delim = FileDelimiters.inferFileDelimiter(inputFile);

            String header = lines.get(0);
            lines.remove(0);
            Map<String,Integer> fieldIndexMap = FileReaderUtils.createFieldsIndexMap(header, delim);

            if(!fieldIndexMap.containsKey(FLD_CHROMOSOME)
            || (!(fieldIndexMap.containsKey(FLD_POSITION_START) && fieldIndexMap.containsKey(FLD_POSITION_END))
                    && !fieldIndexMap.containsKey(FLD_POSITION)))
            {
                GU_LOGGER.error("input file missing required fields: {},{} & {} or {}",
                        FLD_CHROMOSOME, FLD_POSITION_START, FLD_POSITION_END, FLD_POSITION);
                System.exit(1);
            }

            int chrIndex = fieldIndexMap.get(FLD_CHROMOSOME);
            Integer posStartIndex = fieldIndexMap.get(FLD_POSITION_START);
            Integer posEndIndex = fieldIndexMap.get(FLD_POSITION_END);
            Integer posIndex = fieldIndexMap.get(FLD_POSITION);

            BufferedWriter writer = createBufferedWriter(outputFile);

            writer.write(header);
            writer.newLine();

            for(String line : lines)
            {
                String[] values = line.split(delim, -1);
                String chromosome = values[chrIndex];

                StringJoiner sj = new StringJoiner(delim);

                for(int i = 0; i < values.length; ++i)
                {
                    String value = values[i];

                    if(isField(i, posIndex) || isField(i, posStartIndex) || isField(i, posEndIndex))
                    {
                        int origPosition = Integer.parseInt(value);
                        int newPosition = liftoverCache.convertPosition(chromosome, origPosition, destVersion);

                        if(newPosition == UNMAPPED_POSITION)
                        {
                            GU_LOGGER.debug("unmapped location({}:{})", chromosome, origPosition);
                            continue;
                        }

                        sj.add(String.valueOf(newPosition));
                    }
                    else if(i == chrIndex)
                    {
                        sj.add(destVersion.versionedChromosome(chromosome));
                    }
                    else
                    {
                        sj.add(value);
                    }
                }

                writer.write(sj.toString());
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to read/write: {}", e.toString());
            System.exit(1);
        }
    }

    private static boolean isField(int index, @Nullable Integer posIndex)
    {
        return posIndex != null && index == posIndex;
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addPath(INPUT_FILE, true,
                format("Input file (TSV or CSV) - must contain: %s and %s/%s or %s",
                        FLD_CHROMOSOME, FLD_POSITION_START, FLD_POSITION_END, FLD_POSITION));

        configBuilder.addConfigItem(SOURCE_REF_GENOME_VERSION, true, "Source ref genome version");
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output filename");
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        new FileLiftover(configBuilder);
    }
}
