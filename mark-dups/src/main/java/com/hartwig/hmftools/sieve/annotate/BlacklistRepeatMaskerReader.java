package com.hartwig.hmftools.sieve.annotate;

import static com.hartwig.hmftools.sieve.annotate.AnnotateConfig.MD_LOGGER;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import org.jetbrains.annotations.NotNull;

// TODO(m_cooper): Code duplication?
public class BlacklistRepeatMaskerReader
{
    private static final String DELIMITER = ",";

    @NotNull
    static public List<AnnotatedBlacklistRegion> readFromFile(@NotNull final String filepath)
    {
        final Map<BlacklistRegion, List<RepeatMasker>> regionToMasks = new HashMap<>();

        try
        {
            final BufferedReader reader = new BufferedReader(new FileReader(filepath));

            // Drop header.
            reader.readLine();

            String line;
            while((line = reader.readLine()) != null)
            {
                final String[] fields = line.split(DELIMITER);
                if(fields.length != 12)
                {
                    MD_LOGGER.error("Blacklist repeat masker file {} contains a record with {} fields. There should only be 12 fields. ", filepath, fields.length);
                    System.exit(1);
                }

                final BlacklistRepeatMaskerRecord record = parseFields(fields, filepath);
                final BlacklistRegion blacklistRegion = record.getBlacklistRegion();
                final Optional<RepeatMasker> repeatMasker = record.getRepeatMasker();
                if(!regionToMasks.containsKey(blacklistRegion))
                {
                    regionToMasks.put(blacklistRegion, new ArrayList<>());
                }

                if(repeatMasker.isPresent())
                {
                    regionToMasks.get(blacklistRegion).add(repeatMasker.get());
                }
            }

            reader.close();

        }
        catch(Exception e)
        {
            MD_LOGGER.error("An exception was raised while reading the blacklist repeat masker file {}: {}", filepath, e.toString());
            System.exit(1);
        }

        List<AnnotatedBlacklistRegion> output = new ArrayList<>();
        for(var entry : regionToMasks.entrySet())
        {
            final AnnotatedBlacklistRegion annotated = new AnnotatedBlacklistRegion(entry.getKey());
            for(var repeatMasker : entry.getValue())
            {
                annotated.addRepeatMasker(repeatMasker);
            }

            output.add(annotated);
        }

        return output;
    }

    @NotNull
    static private BlacklistRepeatMaskerRecord parseFields(@NotNull final String[] fields, @NotNull final String filepath)
    {
        String chromosome = fields[0];
        int posStart = parseInt(fields[1], "PosStart", filepath);
        int posEnd = parseInt(fields[2], "PosEnd", filepath);
        int sampleCount = parseInt(fields[3], "SampleCount", filepath);
        int depthMin = parseInt(fields[4], "DepthMin", filepath);
        int depthMax = parseInt(fields[5], "DepthMax", filepath);
        Optional<String> repeatType = parseOptionalString(fields[6]);
        Optional<String> repeatInfo = parseOptionalString(fields[7]);
        Optional<Integer> repeatPosStart = parseOptionalInt(fields[8], "RepeatPosStart", filepath);
        Optional<Integer> repeatPosEnd = parseOptionalInt(fields[9], "RepeatPosEnd", filepath);
        Optional<Integer> count = parseOptionalInt(fields[10], "Count", filepath);
        Optional<String> otherInfo = parseOptionalString(fields[11]);

        return new BlacklistRepeatMaskerRecord(chromosome, posStart, posEnd, sampleCount, depthMin, depthMax, repeatType, repeatInfo, repeatPosStart, repeatPosEnd, count, otherInfo);
    }

    static private int parseInt(@NotNull final String str, @NotNull final String fieldName, @NotNull final String filepath)
    {
        try
        {
            return Integer.parseInt(str);
        }
        catch(NumberFormatException e)
        {
            MD_LOGGER.error("While reading a record in the blacklist repeat masker file {}, failed to parse the '{}' field to an int.", filepath, fieldName);
            System.exit(1);
        }

        return 0;
    }

    @NotNull
    static private Optional<Integer> parseOptionalInt(@NotNull final String str, @NotNull final String fieldName,
            @NotNull final String filepath)
    {
        if(str.equals("NA"))
        {
            return Optional.empty();
        }

        try
        {
            return Optional.of(Integer.parseInt(str));
        }
        catch(NumberFormatException e)
        {
            MD_LOGGER.error("While reading a record in the blacklist repeat masker file {}, failed to parse the '{}' field to an int or 'NA'.", filepath, fieldName);
            System.exit(1);
        }

        return Optional.empty();
    }

    @NotNull
    static private Optional<String> parseOptionalString(@NotNull final String str)
    {
        if(str.equals("NA"))
        {
            return Optional.empty();
        }

        return Optional.of(str);
    }
}
