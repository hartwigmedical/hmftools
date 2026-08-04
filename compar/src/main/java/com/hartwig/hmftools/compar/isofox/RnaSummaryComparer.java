package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.rna.RnaStatisticFile;
import com.hartwig.hmftools.common.rna.RnaStatistics;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class RnaSummaryComparer extends ItemComparer
{
    protected enum Fields
    {
        QcStatus,
        TotalFragments,
        DuplicateFragments,
        SplicedFragmentPerc,
        UnsplicedFragmentPerc,
        AltFragmentPerc,
        ChimericFragmentPerc,
        SplicedGeneCount,
        // ReadLength  // a property of the sample, not our tools
        FragLength5th,
        FragLength50th,
        FragLength95th,
        // MedianGCRatio - gone with v2.1 anyway
        MedianGCRatio,
        ForwardStrandPercent;
    }

    public RnaSummaryComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.QcStatus.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.QcStatus.toString()), null));

        mFields.add(new FieldInfo(
                Fields.TotalFragments.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.TotalFragments.toString(), 10.0, 0.01),
                null));

        mFields.add(new FieldInfo(
                Fields.DuplicateFragments.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.DuplicateFragments.toString(), 10.0, 0.01),
                null));

        mFields.add(new FieldInfo(
                Fields.SplicedFragmentPerc.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.SplicedFragmentPerc.toString(), 0.05, null),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.UnsplicedFragmentPerc.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.UnsplicedFragmentPerc.toString(), 0.05, null),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.AltFragmentPerc.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.AltFragmentPerc.toString(), 0.05, null),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.ChimericFragmentPerc.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.ChimericFragmentPerc.toString(), 0.05, null),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.SplicedGeneCount.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.SplicedGeneCount.toString(), 10.0, 0.01),
                null));

        mFields.add(new FieldInfo(
                Fields.FragLength5th.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.FragLength5th.toString(), null, 0.05),
                "%.1f"));

        mFields.add(new FieldInfo(
                Fields.FragLength50th.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.FragLength50th.toString(), null, 0.05),
                "%.1f"));

        mFields.add(new FieldInfo(
                Fields.FragLength95th.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.FragLength95th.toString(), null, 0.05),
                "%.1f"));

        mFields.add(new FieldInfo(
                Fields.MedianGCRatio.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.MedianGCRatio.toString(), 0.01, null),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.ForwardStrandPercent.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.ForwardStrandPercent.toString(), 0.01, null),
                "%.2f"));
    }

    @Override
    public CategoryType category()
    {
        return CategoryType.RNA_SUMMARY;
    }

    @Override
    public boolean hasReportable()
    {
        return false;
    }

    @Override
    public List<String> displayFieldNames()
    {
        return Arrays.stream(RnaSummaryComparer.Fields.values()).map(x -> x.toString()).collect(Collectors.toList());
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(determineFileName(sampleId, fileSources)));
            RnaStatistics rnaStatistics = RnaStatisticFile.fromLines(lines);
            comparableItems.add(new RnaSummaryData(rnaStatistics, mFields));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Isofox Summary data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }

    private static String determineFileName(final String sampleId, final PipelineSourcePaths fileSources)
    {
        String filename = RnaStatisticFile.generateFilename(fileSources.Isofox, sampleId);
        String oldFilename = filename.replace(TSV_EXTENSION, CSV_EXTENSION);

        if(!Files.exists(Paths.get(filename)) && Files.exists(Paths.get(oldFilename)))
        {
            return oldFilename;
        }
        else
        {
            return filename;
        }
    }
}
