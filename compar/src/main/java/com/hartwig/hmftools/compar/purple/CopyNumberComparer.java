package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.common.CategoryType.COPY_NUMBER;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.FieldCheckCache.getOrMakeFieldCheck;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class CopyNumberComparer extends ItemComparer
{
    protected enum Fields
    {
        CopyNumber,
        MajorAlleleCopyNumber,
        Method;
    }

    public CopyNumberComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.CopyNumber.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.CopyNumber.toString(), 0.5, 0.15),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.MajorAlleleCopyNumber.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.MajorAlleleCopyNumber.toString(), 0.5, 0.15),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.Method.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.Method.toString()), null));
    }

    @Override
    public CategoryType category() { return COPY_NUMBER; }

    @Override
    public boolean hasReportable() { return false; }

    @Override
    public List<String> displayFieldNames()
    {
        return Arrays.stream(Fields.values()).map(x -> x.toString()).collect(Collectors.toList());
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            List<PurpleCopyNumber> copyNumbers = PurpleCopyNumberFile.read(PurpleCopyNumberFile.generateFilenameForReading(
                    fileSources.Purple, sampleId));

            copyNumbers.forEach(x -> comparableItems.add(createCopyNumberData(x)));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to read copy number data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }

    private CopyNumberData createCopyNumberData(final PurpleCopyNumber copyNumber)
    {
        return new CopyNumberData(
                copyNumber.chromosome(), copyNumber.start(), copyNumber.end(),
                copyNumber.averageTumorCopyNumber(), copyNumber.majorAlleleCopyNumber(),
                copyNumber.method(), mFields);
    }
}
