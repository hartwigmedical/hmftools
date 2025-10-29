package com.hartwig.hmftools.panelbuilder;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public record CustomSv(
        BasePosition startPosition,
        Orientation startOrientation,
        BasePosition endPosition,
        Orientation endOrientation,
        String insertSequence,
        // Arbitrary descriptor for the user.
        String extraInfo
)
{
    private enum FileFields
    {
        ChromosomeStart,
        PositionStart,
        OrientationStart,
        ChromosomeEnd,
        PositionEnd,
        OrientationEnd,
        InsertSequence,
        ExtraInfo;
    }

    private static final Logger LOGGER = LogManager.getLogger(CustomSv.class);

    public static List<CustomSv> readFromFile(final String filePath)
    {
        LOGGER.debug("Reading custom structural variants from file: {}", filePath);

        try(DelimFileReader reader = new DelimFileReader(filePath))
        {
            List<CustomSv> customSvs = reader.stream().map(row ->
            {
                BasePosition startPosition = new BasePosition(row.get(FileFields.ChromosomeStart), row.getInt(FileFields.PositionStart));
                Orientation startOrientation = Orientation.fromByteStr(row.get(FileFields.OrientationStart));
                BasePosition endPosition = new BasePosition(row.get(FileFields.ChromosomeEnd), row.getInt(FileFields.PositionEnd));
                Orientation endOrientation = Orientation.fromByteStr(row.get(FileFields.OrientationEnd));
                String insertSequence = row.getOrNull(FileFields.InsertSequence);
                if(insertSequence == null)
                {
                    insertSequence = "";
                }
                String extraInfo = row.get(FileFields.ExtraInfo);
                return new CustomSv(startPosition, startOrientation, endPosition, endOrientation, insertSequence, extraInfo);
            }).toList();

            LOGGER.debug("Read {} custom structural variants from {}", customSvs.size(), filePath);
            return customSvs;
        }
    }

    public static void writeToFile(final List<CustomSv> customSvs, final String filePath)
    {
        LOGGER.debug("Writing custom structural variants to file: {}", filePath);

        try(DelimFileWriter<CustomSv> writer = new DelimFileWriter<>(filePath, FileFields.values(), CustomSv::writeObj))
        {
            customSvs.forEach(writer::writeRow);
        }
    }

    private static void writeObj(final CustomSv sv, final DelimFileWriter.Row row)
    {
        row.set(FileFields.ChromosomeStart, sv.startPosition().Chromosome);
        row.set(FileFields.PositionStart, sv.startPosition().Position);
        row.set(FileFields.OrientationStart, Byte.toString(sv.startOrientation().asByte()));
        row.set(FileFields.ChromosomeEnd, sv.endPosition().Chromosome);
        row.set(FileFields.PositionEnd, sv.endPosition().Position);
        row.set(FileFields.OrientationEnd, Byte.toString(sv.endOrientation().asByte()));
        row.set(FileFields.InsertSequence, sv.insertSequence());
        row.set(FileFields.ExtraInfo, sv.extraInfo());
    }
}
