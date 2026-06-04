package com.hartwig.hmftools.panelbuilder;

import java.util.List;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

// User-specified small insertion/deletion variant.
// Represented as a ref sequence which is replaced by an alt sequence, at a specific position.
public record CustomSmallVariant(
        BasePosition position,
        // Ref sequence is assumed to match the ref genome at the position. Technically, we only need the length of this sequence.
        // We specify it like this because it's what a typical user would expect (similar to VCF format).
        String refSequence,
        String altSequence,
        // Arbitrary descriptor for the user.
        String extraInfo,
        // If null, use the default quality score minimum.
        @Nullable Double qualityScoreMin
)
{
    public CustomSmallVariant
    {
        if(refSequence.isEmpty())
        {
            // The first base of refSequence must match the position, so empty is not allowed.
            // If you want to insert 1 base, need to set refSequence=A and altSequence=AB
            throw new IllegalArgumentException("refSequence cannot be empty");
        }
        if(refSequence.equals(altSequence))
        {
            throw new IllegalArgumentException("refSequence and altSequence imply no change");
        }
    }

    private enum Columns
    {
        Chromosome,
        Position,
        RefSequence,
        AltSequence,
        ExtraInfo,
        QualityScoreMin
    }

    private static final Logger LOGGER = LogManager.getLogger(CustomSmallVariant.class);

    public static List<CustomSmallVariant> readFromFile(final String filePath)
    {
        LOGGER.debug("Reading custom small variants from file: {}", filePath);

        try(DelimFileReader reader = new DelimFileReader(filePath))
        {
            List<CustomSmallVariant> items = reader.stream().map(row ->
            {
                BasePosition position = new BasePosition(row.getString(Columns.Chromosome), row.getInt(Columns.Position));
                String refSequence = row.getString(Columns.RefSequence);
                String altSequence = row.getString(Columns.AltSequence);
                String extraInfo = row.getString(Columns.ExtraInfo);
                Double qualityScoreMin = row.getDoubleOrNull(Columns.QualityScoreMin);
                return new CustomSmallVariant(position, refSequence, altSequence, extraInfo, qualityScoreMin);
            }).toList();

            LOGGER.debug("Read {} custom small variants from {}", items.size(), filePath);
            return items;
        }
    }

    public static void writeToFile(final List<CustomSmallVariant> items, final String filePath)
    {
        LOGGER.debug("Writing custom small variants to file: {}", filePath);

        try(DelimFileWriter<CustomSmallVariant> writer = new DelimFileWriter<>(filePath, Columns.values(), CustomSmallVariant::writeObj))
        {
            items.forEach(writer::writeRow);
        }
    }

    private static void writeObj(final CustomSmallVariant obj, final DelimFileWriter.Row row)
    {
        row.set(Columns.Chromosome, obj.position().Chromosome);
        row.set(Columns.Position, obj.position().Position);
        row.set(Columns.RefSequence, obj.refSequence());
        row.set(Columns.AltSequence, obj.altSequence());
        row.set(Columns.ExtraInfo, obj.extraInfo());
        row.setOrNull(Columns.QualityScoreMin, obj.qualityScoreMin());
    }
}
