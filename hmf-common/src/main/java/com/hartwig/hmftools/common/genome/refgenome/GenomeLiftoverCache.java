package com.hartwig.hmftools.common.genome.refgenome;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.LOGGER;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.stripChrPrefix;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class GenomeLiftoverCache
{
    private final Map<String,List<CoordMapping>> mMappings;

    public static final String LIFTOVER_MAPPING_FILE = "liftover_mapping";
    public static final String LIFTOVER_MAPPING_FILE_DESC = "Liftover mapping file";

    public static final int UNMAPPED_POSITION = -1;

    private static final String NEG_ORIENT = "-";

    public GenomeLiftoverCache() { this(false); }

    public GenomeLiftoverCache(boolean loadDefaultMappings)
    {
        mMappings = Maps.newHashMap();

        if(loadDefaultMappings)
            loadDefaultResource();
    }

    public boolean hasMappings() { return !mMappings.isEmpty(); }

    public List<CoordMapping> getChromosomeMappings(final String chromosome) { return mMappings.get(chromosome); }

    public int convertPosition(final String chromosome, final int position) { return convertPositionTo38(chromosome, position); }

    public int convertPosition(final String chromosome, final int position, final RefGenomeVersion destinationVersion)
    {
        return destinationVersion == V38 ? convertPositionTo38(chromosome, position) : convertPositionTo37(chromosome, position);
    }

    public int convertPositionTo37(final String chromosome, final int position)
    {
        List<CoordMapping> mappings = mMappings.get(stripChrPrefix(chromosome));

        if(mappings == null || mappings.isEmpty())
            return position;

        for(CoordMapping mapping : mappings)
        {
            if(mapping.DestEnd < position)
                continue;

            if(mapping.DestStart > position)
                return UNMAPPED_POSITION;

            return mapping.reversePosition(position);
        }

        return UNMAPPED_POSITION;
    }

    public int convertPositionTo38(final String chromosome, final int position)
    {
        List<CoordMapping> mappings = mMappings.get(chromosome);

        if(mappings == null || mappings.isEmpty())
            return position;

        for(CoordMapping mapping : mappings)
        {
            if(mapping.SourceEnd < position)
                continue;

            if(mapping.SourceStart > position)
                return UNMAPPED_POSITION;

            return mapping.convertPosition(position);
        }

        return UNMAPPED_POSITION;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(LIFTOVER_MAPPING_FILE, true, LIFTOVER_MAPPING_FILE_DESC);
    }

    public boolean loadFile(final String filename)
    {
        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            if(loadFile(lines))
            {
                LOGGER.info("loaded {} genome liftover mapping entries from file: {}", mMappings.size(), filename);
                return true;
            }
        }
        catch(Exception e)
        {
            LOGGER.error("failed to load genome liftover entries from file: {}", filename, e.toString());
        }

        return false;
    }

    private void loadDefaultResource()
    {
        final InputStream inputStream = GenomeLiftoverCache.class.getResourceAsStream("/refgenome/hg37_38_mapping.tsv");
        List<String> lines = new BufferedReader(new InputStreamReader(inputStream)).lines().collect(Collectors.toList());
        loadFile(lines);
    }

    public boolean loadFile(final List<String> lines)
    {
        String header = lines.get(0);
        lines.remove(0);

        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

        List<CoordMapping> chrMappings = null;
        String currentChromosome = "";

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            if(values.length < 8)
            {
                LOGGER.error("invalid mapping line: {}", line);
                return false;
            }

            String chromosome = values[0];

            if(!currentChromosome.equals(chromosome))
            {
                currentChromosome = chromosome;
                chrMappings = Lists.newArrayList();
                mMappings.put(chromosome, chrMappings);
            }

            // note +1 for start positions since source file is in BED style
            // Chr37	Start37	End37	Orient37	Chr38 	Start38	End38	Orient38
            boolean isReverse = values[fieldsIndexMap.get("Orient38")].equals(NEG_ORIENT);

            chrMappings.add(new CoordMapping(
                    values[fieldsIndexMap.get("Chr37")],
                    Integer.parseInt(values[fieldsIndexMap.get("Start37")]) + 1,
                    Integer.parseInt(values[fieldsIndexMap.get("End37")]),
                    Integer.parseInt(values[fieldsIndexMap.get("Start38")]) + 1,
                    Integer.parseInt(values[fieldsIndexMap.get("End38")]),
                    isReverse));
        }

        return true;
    }

    @VisibleForTesting
    public void addMapping(final String chromosome, int sourceStart, int sourceEnd, int destStart, int destEnd, boolean reverse)
    {
        List<CoordMapping> chrMappings = mMappings.get(chromosome);

        if(chrMappings == null)
        {
            chrMappings = Lists.newArrayList();
            mMappings.put(chromosome, chrMappings);
        }
        chrMappings.add(new CoordMapping(chromosome, sourceStart, sourceEnd, destStart, destEnd, reverse));
    }
}
