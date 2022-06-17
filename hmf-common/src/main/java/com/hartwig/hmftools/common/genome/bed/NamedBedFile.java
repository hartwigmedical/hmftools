package com.hartwig.hmftools.common.genome.bed;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;
import java.util.zip.GZIPOutputStream;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.readers.LineIterator;

public final class NamedBedFile
{
    private static final Logger LOGGER = LogManager.getLogger(NamedBedFile.class);
    private static final String DELIMITER = "\t";

    private NamedBedFile()
    {
    }

    public static void writeUnnamedBedFile(final String filename, final List<GenomeRegion> regions) throws IOException
    {
        List<String> strings = regions.stream().map(NamedBedFile::asBed).collect(Collectors.toList());
        write(filename, strings);
    }

    public static void writeBedFile(final String filename, final List<NamedBed> regions) throws IOException
    {
        List<String> strings = regions.stream().map(NamedBedFile::asBed).collect(Collectors.toList());
        write(filename, strings);
    }

    public static void write(final String filename, final List<String> lines) throws IOException
    {
        try(FileOutputStream output = new FileOutputStream(filename))
        {
            OutputStream transformedOutput = filename.endsWith(".gz") ? new GZIPOutputStream(output) : output;
            try(Writer writer = new OutputStreamWriter(transformedOutput, StandardCharsets.UTF_8))
            {
                for(CharSequence line : lines)
                {
                    writer.append(line);
                    writer.append("\n");
                }
            }
        }
    }

    @NotNull
    public static List<NamedBed> readBedFile(String bedFile) throws IOException
    {
        List<NamedBed> result = Lists.newArrayList();
        NamedBed prevRegion = null;
        try(final AbstractFeatureReader<BEDFeature, LineIterator> reader = getFeatureReader(bedFile, new BEDCodec(), false))
        {
            for(final BEDFeature bedFeature : reader.iterator())
            {
                final NamedBed namedBed = fromBedFeature(bedFeature);
                if(namedBed.end() < namedBed.start())
                {
                    LOGGER.warn("Invalid genome region found in chromosome {}: start={}, end={}",
                            namedBed.chromosome(),
                            namedBed.start(),
                            namedBed.end());
                }
                else
                {
                    if(prevRegion != null && namedBed.chromosome().equals(prevRegion.chromosome())
                            && prevRegion.end() >= namedBed.start())
                    {
                        LOGGER.warn("BED file is not sorted, please fix! Current={}, Previous={}", namedBed, prevRegion);
                    }
                    else
                    {
                        result.add(namedBed);
                        prevRegion = namedBed;
                    }
                }
            }
        }

        return result;
    }

    static NamedBed fromBedFeature(BEDFeature feature)
    {
        String name = feature.getName();
        return ImmutableNamedBed.builder()
                .chromosome(feature.getContig())
                .start(feature.getStart())
                .end(feature.getEnd())
                .name(name == null ? Strings.EMPTY : name)
                .build();
    }

    private static String asBed(final GenomeRegion region)
    {
        return new StringJoiner(DELIMITER).add(region.chromosome())
                .add(String.valueOf(region.start() - 1))
                .add(String.valueOf(region.end()))
                .toString();
    }

    private static String asBed(final NamedBed region)
    {
        return new StringJoiner(DELIMITER).add(region.chromosome())
                .add(String.valueOf(region.start() - 1))
                .add(String.valueOf(region.end()))
                .add(region.name())
                .toString();
    }
}
