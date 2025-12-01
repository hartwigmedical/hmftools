package com.hartwig.hmftools.panelbuilder;

import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

public record CustomRegion(
        ChrBaseRegion region,
        // Arbitrary descriptor for the user.
        String extraInfo,
        // If null, use the default quality score minimum.
        @Nullable Double qualityScoreMin
)
{
    private enum Columns
    {
        Chromosome,
        PositionStart,
        PositionEnd,
        ExtraInfo,
        QualityScoreMin
    }

    private static final Logger LOGGER = LogManager.getLogger(CustomRegion.class);

    public static List<CustomRegion> readFromFile(final String filePath)
    {
        LOGGER.debug("Reading custom regions from file: {}", filePath);

        try(DelimFileReader reader = new DelimFileReader(filePath))
        {
            List<CustomRegion> customRegions = reader.stream().map(row ->
            {
                String chromosome = row.get(Columns.Chromosome);
                int startPosition = row.getInt(Columns.PositionStart);
                int endPosition = row.getInt(Columns.PositionEnd);
                String extraInfo = row.get(Columns.ExtraInfo);
                Double qualityScoreMin = row.getDoubleOrNull(Columns.QualityScoreMin);
                return new CustomRegion(new ChrBaseRegion(chromosome, startPosition, endPosition), extraInfo, qualityScoreMin);
            }).toList();

            LOGGER.debug("Read {} custom regions from {}", customRegions.size(), filePath);
            return customRegions;
        }
    }

    public static void writeToFile(final List<CustomRegion> customRegions, final String filePath)
    {
        LOGGER.debug("Writing custom regions to file: {}", filePath);

        try(DelimFileWriter<CustomRegion> writer = new DelimFileWriter<>(filePath, Columns.values(), CustomRegion::writeObj))
        {
            customRegions.forEach(writer::writeRow);
        }
    }

    private static void writeObj(final CustomRegion customRegion, final DelimFileWriter.Row row)
    {
        row.set(Columns.Chromosome, customRegion.region().chromosome());
        row.set(Columns.PositionStart, customRegion.region().start());
        row.set(Columns.PositionEnd, customRegion.region().end());
        row.set(Columns.ExtraInfo, customRegion.extraInfo());
        row.setOrNull(Columns.QualityScoreMin, customRegion.qualityScoreMin());
    }
}
