package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.common.CategoryType.GENE_COPY_NUMBER;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class GeneCopyNumberComparer extends ItemComparer
{
    protected enum Fields
    {
        MinCopyNumber,
        MaxCopyNumber,
        MinRegionStart,
        MinRegionEnd;
    }

    public GeneCopyNumberComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.MinCopyNumber.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.MinCopyNumber.toString(), 0.5, 0.15),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.MaxCopyNumber.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.MaxCopyNumber.toString(), 0.5, 0.15),
                "%.2f"));

        mFields.add(FieldInfo.displayOnly(Fields.MinRegionStart.toString(), null));
        mFields.add(FieldInfo.displayOnly(Fields.MinRegionEnd.toString(), null));
    }

    @Override
    public CategoryType category() { return GENE_COPY_NUMBER; }

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
        final List<ComparableItem> items = Lists.newArrayList();

        // load only the genes in the driver catalog if the driver category is being run

        final Set<String> driverGenes = mConfig.DriverGenes;

        if(driverGenes.isEmpty())
            return items;

        try
        {
            List<GeneCopyNumber> copyNumbers = GeneCopyNumberFile.read(GeneCopyNumberFile.generateFilename(
                    fileSources.Purple, sampleId));

            copyNumbers.stream().filter(x -> driverGenes.contains(x.geneName())).forEach(x -> items.add(new GeneCopyNumberData(x, mFields)));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to read gene copy number data: {}", sampleId, e.toString());
            return null;
        }

        return items;
    }
}
