package com.hartwig.hmftools.panelbuilder;

import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;

import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public record CustomRegion(
        ChrBaseRegion region,
        // Arbitrary descriptor for the user.
        String extraInfo
)
{
    private static final String FLD_EXTRA_INFO = "ExtraInfo";

    private static final List<String> COLUMNS = List.of(FLD_CHROMOSOME, FLD_POSITION_START, FLD_POSITION_END, FLD_EXTRA_INFO);

    private static final Logger LOGGER = LogManager.getLogger(CustomRegion.class);

    public static List<CustomRegion> readFromFile(final String filePath)
    {
        LOGGER.debug("Reading custom regions from file: {}", filePath);

        try(DelimFileReader reader = new DelimFileReader(filePath))
        {
            int chromosomeIdx = requireNonNull(reader.getColumnIndex(FLD_CHROMOSOME));
            int posStartIdx = requireNonNull(reader.getColumnIndex(FLD_POSITION_START));
            int posEndIdx = requireNonNull(reader.getColumnIndex(FLD_POSITION_END));
            int extraInfoIdx = requireNonNull(reader.getColumnIndex(FLD_EXTRA_INFO));

            List<CustomRegion> regions = reader.stream().map(row ->
            {
                String chromosome = row.get(chromosomeIdx);
                int start = row.getInt(posStartIdx);
                int end = row.getInt(posEndIdx);
                String extraInfo = row.get(extraInfoIdx);
                ChrBaseRegion baseRegion = new ChrBaseRegion(chromosome, start, end);
                return new CustomRegion(baseRegion, extraInfo);
            }).toList();

            LOGGER.debug("Read {} custom regions from {}", regions.size(), filePath);
            return regions;
        }
    }

    public static void writeToFile(final List<CustomRegion> regions, final String filePath)
    {
        LOGGER.debug("Writing custom regions to file: {}", filePath);

        try(DelimFileWriter<CustomRegion> writer = new DelimFileWriter<>(filePath, COLUMNS, CustomRegion::writeObj))
        {
            regions.forEach(writer::writeRow);
        }
    }

    private static void writeObj(final CustomRegion region, final DelimFileWriter.Row row)
    {
        row.set(FLD_CHROMOSOME, region.region().chromosome());
        row.set(FLD_POSITION_START, region.region().start());
        row.set(FLD_POSITION_END, region.region().end());
        row.set(FLD_EXTRA_INFO, region.extraInfo());
    }
}
