package com.hartwig.hmftools.geneutils.utils;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
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

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.file.FileDelimiters;
import com.hartwig.hmftools.common.utils.file.FileReaderUtils;

import org.jetbrains.annotations.NotNull;

public class GenerateProbes
{
    private static final String INPUT_FILE = "input_file";
    private static final String OUTPUT_FILE = "output_file";
    private static final String PROBE_LENGTH = "probe_length";

    public GenerateProbes(final ConfigBuilder configBuilder)
    {
        String inputFile = configBuilder.getValue(INPUT_FILE);
        String outputFile = configBuilder.getValue(OUTPUT_FILE);
        int probeLength = configBuilder.getInteger(PROBE_LENGTH);

        RefGenomeInterface refGenome = loadRefGenome(configBuilder.getValue(REF_GENOME));

        GU_LOGGER.info("generating probes length({}) for file({}) to output({})",
                probeLength, inputFile, outputFile);

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(inputFile));
            String delim = FileDelimiters.inferFileDelimiter(inputFile);

            String header = lines.get(0);
            lines.remove(0);
            Map<String,Integer> fieldIndexMap = FileReaderUtils.createFieldsIndexMap(header, delim);

            if(!fieldIndexMap.containsKey(FLD_CHROMOSOME) || !fieldIndexMap.containsKey(FLD_POSITION))
            {
                GU_LOGGER.error("input file missing required fields: {}, {}", FLD_CHROMOSOME, FLD_POSITION);
                System.exit(1);
            }

            int chrIndex = fieldIndexMap.get(FLD_CHROMOSOME);
            int posIndex = fieldIndexMap.get(FLD_POSITION);

            BufferedWriter writer = createBufferedWriter(outputFile);

            StringJoiner sj = new StringJoiner(delim);
            sj.add(header);
            sj.add("Probe");
            sj.add("ProbeStart");
            sj.add("ProbeEnd");
            writer.write(sj.toString());
            writer.newLine();

            int probeLengthLeft = probeLength / 2;
            int probeLengthRight = probeLength  - probeLengthLeft;

            for(String line : lines)
            {
                String[] values = line.split(delim, -1);
                String chromosome = values[chrIndex];
                int position = Integer.parseInt(values[posIndex]);

                int probeStart = position - probeLengthLeft;
                int probeEnd = position + probeLengthRight - 1;

                String probeSequence = refGenome.getBaseString(chromosome, probeStart, probeEnd);

                sj = new StringJoiner(delim);
                sj.add(line);
                sj.add(probeSequence);
                sj.add(String.valueOf(probeStart));
                sj.add(String.valueOf(probeEnd));
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

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addPath(INPUT_FILE, true, "Input file with probe locations");
        configBuilder.addPath(REF_GENOME, true, REF_GENOME_CFG_DESC);
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output filename");
        configBuilder.addRequiredInteger(PROBE_LENGTH, "Probe length");
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        new GenerateProbes(configBuilder);
    }
}
