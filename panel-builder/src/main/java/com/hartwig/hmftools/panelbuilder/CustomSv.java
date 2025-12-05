package com.hartwig.hmftools.panelbuilder;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

public record CustomSv(
        BasePosition startPosition,
        Orientation startOrientation,
        BasePosition endPosition,
        Orientation endOrientation,
        String insertSequence,
        // Arbitrary descriptor for the user.
        String extraInfo,
        // If null, use the default quality score minimum.
        @Nullable Double qualityScoreMin
)
{
    private enum Columns
    {
        ChromosomeStart,
        PositionStart,
        OrientationStart,
        ChromosomeEnd,
        PositionEnd,
        OrientationEnd,
        InsertSequence,
        ExtraInfo,
        QualityScoreMin
    }

    private static final Logger LOGGER = LogManager.getLogger(CustomSv.class);

    public static List<CustomSv> readFromFile(final String filePath)
    {
        LOGGER.debug("Reading custom structural variants from file: {}", filePath);

        try(DelimFileReader reader = new DelimFileReader(filePath))
        {
            List<CustomSv> customSvs = reader.stream().map(row ->
            {
                BasePosition startPosition = new BasePosition(row.getString(Columns.ChromosomeStart), row.getInt(Columns.PositionStart));
                Orientation startOrientation = Orientation.fromByteStr(row.getString(Columns.OrientationStart));
                BasePosition endPosition = new BasePosition(row.getString(Columns.ChromosomeEnd), row.getInt(Columns.PositionEnd));
                Orientation endOrientation = Orientation.fromByteStr(row.getString(Columns.OrientationEnd));
                String insertSequence = row.getStringOrNull(Columns.InsertSequence);
                if(insertSequence == null)
                {
                    insertSequence = "";
                }
                String extraInfo = row.getString(Columns.ExtraInfo);
                Double qualityScoreMin = row.getDoubleOrNull(Columns.QualityScoreMin);
                return new CustomSv(startPosition, startOrientation, endPosition, endOrientation, insertSequence, extraInfo, qualityScoreMin);
            }).toList();

            LOGGER.debug("Read {} custom structural variants from {}", customSvs.size(), filePath);
            return customSvs;
        }
    }

    public static void writeToFile(final List<CustomSv> customSvs, final String filePath)
    {
        LOGGER.debug("Writing custom structural variants to file: {}", filePath);

        try(DelimFileWriter<CustomSv> writer = new DelimFileWriter<>(filePath, Columns.values(), CustomSv::writeObj))
        {
            customSvs.forEach(writer::writeRow);
        }
    }

    private static void writeObj(final CustomSv sv, final DelimFileWriter.Row row)
    {
        row.set(Columns.ChromosomeStart, sv.startPosition().Chromosome);
        row.set(Columns.PositionStart, sv.startPosition().Position);
        row.set(Columns.OrientationStart, Byte.toString(sv.startOrientation().asByte()));
        row.set(Columns.ChromosomeEnd, sv.endPosition().Chromosome);
        row.set(Columns.PositionEnd, sv.endPosition().Position);
        row.set(Columns.OrientationEnd, Byte.toString(sv.endOrientation().asByte()));
        row.set(Columns.InsertSequence, sv.insertSequence());
        row.set(Columns.ExtraInfo, sv.extraInfo());
        row.setOrNull(Columns.QualityScoreMin, sv.qualityScoreMin());
    }
}
