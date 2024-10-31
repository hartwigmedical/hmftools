package com.hartwig.hmftools.geneutils.utils;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache.UNMAPPED_POSITION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.INVALID_FIELD;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.getChromosomeFieldIndex;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.getPositionEndFieldIndex;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.getPositionFieldIndex;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.getPositionStartFieldIndex;
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
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.file.FileDelimiters;
import com.hartwig.hmftools.common.utils.file.FileReaderUtils;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class FileLiftover
{
    private static final String INPUT_FILE = "input_file";
    private static final String OUTPUT_FILE = "output_file";
    private static final String SORT_OUTPUT = "sort_output";
    private static final String SOURCE_REF_GENOME_VERSION = "source_ref_genome_version";

    public FileLiftover(final ConfigBuilder configBuilder)
    {
        String inputFile = configBuilder.getValue(INPUT_FILE);
        String outputFile = configBuilder.getValue(OUTPUT_FILE);

        RefGenomeVersion sourceVersion = RefGenomeVersion.from(configBuilder.getValue(SOURCE_REF_GENOME_VERSION));
        RefGenomeVersion destVersion = sourceVersion.is37() ? V38 : V37;

        boolean sortOutput = configBuilder.hasFlag(SORT_OUTPUT);

        GenomeLiftoverCache liftoverCache = new GenomeLiftoverCache(true);

        GU_LOGGER.info("lifting over file({}) to ref genome version({})", inputFile, destVersion);

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(inputFile));
            String delim = FileDelimiters.inferFileDelimiter(inputFile);

            String header = lines.get(0);
            lines.remove(0);
            Map<String,Integer> fieldIndexMap = FileReaderUtils.createFieldsIndexMap(header, delim);

            int chrIndex = getChromosomeFieldIndex(fieldIndexMap);
            int posIndex = getPositionFieldIndex(fieldIndexMap);
            Integer posStartIndex = getPositionStartFieldIndex(fieldIndexMap);
            Integer posEndIndex = getPositionEndFieldIndex(fieldIndexMap);

            boolean fieldsPresent = chrIndex != INVALID_FIELD
                    && (posIndex != INVALID_FIELD || (posStartIndex != INVALID_FIELD && posEndIndex != INVALID_FIELD));

            if(!fieldsPresent)
            {
                GU_LOGGER.error("input file missing required fields: {},{} & {} or {}",
                        FLD_CHROMOSOME, FLD_POSITION_START, FLD_POSITION_END, FLD_POSITION);
                System.exit(1);
            }

            Map<String,List<RegionEntry>> chrRegionEntries = Maps.newHashMap();
            List<String> chromosomes = Lists.newArrayList();

            BufferedWriter writer = createBufferedWriter(outputFile);

            writer.write(header);
            writer.newLine();

            int failedLiftoverCount = 0;
            String currentChromosome = "";
            List<RegionEntry> currentRegions = null;

            for(String line : lines)
            {
                String[] values = line.split(delim, -1);

                String chromosome = values[chrIndex];

                if(!currentChromosome.equals(chromosome))
                {
                    currentChromosome = chromosome;
                    chromosomes.add(chromosome);
                    currentRegions = Lists.newArrayList();
                    chrRegionEntries.put(chromosome, currentRegions);
                }

                boolean liftoverOk = true;
                int origPosStart = 0;
                int origPosEnd = 0;
                int newPosStart = 0;
                int newPosEnd = 0;

                // first establish the coordinates
                for(int i = 0; i < values.length; ++i)
                {
                    String value = values[i];

                    if(isField(i, posIndex) || isField(i, posStartIndex) || isField(i, posEndIndex))
                    {
                        int origPosition = Integer.parseInt(value);
                        int newPosition = liftoverCache.convertPosition(chromosome, origPosition, destVersion);

                        if(newPosition == UNMAPPED_POSITION)
                        {
                            GU_LOGGER.debug("skipped writing unmapped location({}:{})", chromosome, origPosition);
                            liftoverOk = false;
                            break;
                        }

                        if(isField(i, posEndIndex))
                        {
                            newPosEnd = newPosition;
                            origPosEnd = origPosition;
                        }
                        else
                        {
                            newPosStart = newPosition;
                            origPosStart = origPosition;
                        }
                    }
                }

                if(!liftoverOk)
                {
                    ++failedLiftoverCount;
                    continue;
                }

                // check for position reversal if using start and end
                if(newPosEnd > 0 && newPosStart > newPosEnd)
                {
                    GU_LOGGER.info("location src({}:{}-{}) dest({}-{}) has positions reversed",
                            chromosome, origPosStart, origPosEnd, newPosStart, newPosEnd);

                    int tmp = newPosEnd;
                    newPosEnd = newPosStart;
                    newPosStart = tmp;
                }

                StringJoiner sj = new StringJoiner(delim);

                for(int i = 0; i < values.length; ++i)
                {
                    String value = values[i];

                    if(isField(i, posIndex) || isField(i, posStartIndex))
                    {
                        sj.add(String.valueOf(newPosStart));
                    }
                    else if(isField(i, posEndIndex))
                    {
                        sj.add(String.valueOf(newPosEnd));
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

                if(liftoverOk)
                {
                    if(sortOutput)
                    {
                        currentRegions.add(new RegionEntry(newPosStart, newPosEnd, sj.toString()));
                    }
                    else
                    {
                        writer.write(sj.toString());
                        writer.newLine();
                    }
                }
                else
                {
                    ++failedLiftoverCount;
                }
            }

            if(sortOutput)
            {
                for(String chromosome : chromosomes)
                {
                    List<RegionEntry> regions = chrRegionEntries.get(chromosome);

                    Collections.sort(regions);

                    for(RegionEntry region : regions)
                    {
                        writer.write(region.Data);
                        writer.newLine();
                    }
                }
            }

            writer.close();

            if(failedLiftoverCount > 0)
            {
                GU_LOGGER.info("skipped lifting over {} entries", failedLiftoverCount);
            }
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to read/write: {}", e.toString());
            System.exit(1);
        }
    }

    private class RegionEntry extends BaseRegion
    {
        public final String Data;

        public RegionEntry(final int posStart, final int posEnd, final String data)
        {
            super(posStart, posEnd);
            Data = data;
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
        configBuilder.addFlag(SORT_OUTPUT, "Sort output by chromosome and position");
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        new FileLiftover(configBuilder);
    }
}
