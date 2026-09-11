package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.common.CategoryType.GERMLINE_AMP_DEL;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;
import static com.hartwig.hmftools.compar.common.ComparConstants.COPY_NUMBER_ABS_THRESHOLD;
import static com.hartwig.hmftools.compar.common.ComparConstants.COPY_NUMBER_PERC_THRESHOLD;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class GermlineAmpDelComparer extends ItemComparer
{
    protected enum Fields
    {
        Reported,
        GermlineStatus,
        TumorStatus,
        GermlineCopyNumber,
        TumorCopyNumber;
    }

    public GermlineAmpDelComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.Reported.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.Reported.toString()), null));

        mFields.add(new FieldInfo(
                Fields.GermlineStatus.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.GermlineStatus.toString()), null));

        mFields.add(new FieldInfo(
                Fields.TumorStatus.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.TumorStatus.toString()), null));

        mFields.add(new FieldInfo(
                Fields.GermlineCopyNumber.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.GermlineCopyNumber.toString(), COPY_NUMBER_ABS_THRESHOLD, COPY_NUMBER_PERC_THRESHOLD),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.TumorCopyNumber.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.TumorCopyNumber.toString(), COPY_NUMBER_ABS_THRESHOLD, COPY_NUMBER_PERC_THRESHOLD),
                "%.2f"));
    }

    @Override
    public CategoryType category() { return GERMLINE_AMP_DEL; }

    /*
    @Override
    public List<Field> fields(final MatchLevel matchLevel)
    {
        return List.of(
                new StringField(FLD_REPORTED, i -> ((GermlineAmpDelData) i).AmpDelData.Reported.toString(), true),
                new StringField(FLD_GERMLINE_STATUS, i -> ((GermlineAmpDelData) i).AmpDelData.NormalStatus.toString(),
                        true),
                new StringField(FLD_TUMOR_STATUS, i -> ((GermlineAmpDelData) i).AmpDelData.TumorStatus.toString(),
                        true),
                new DoubleField(FLD_GERMLINE_CN, i -> ((GermlineAmpDelData) i).AmpDelData.GermlineCopyNumber,
                        true, 0.2, 0.1, "%.2f"),
                new DoubleField(FLD_TUMOR_CN, i -> ((GermlineAmpDelData) i).AmpDelData.TumorCopyNumber,
                        true, 0.2, 0.1, "%.2f"),
                new StringField(FLD_CHROMOSOME, i -> ((GermlineAmpDelData) i).mComparisonChromosome,
                        true),
                new StringField(FLD_CHROMOSOME_BAND, i -> ((GermlineAmpDelData) i).AmpDelData.ChromosomeBand, true)
        );
    }
    */

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
            String germlineAmpDelFile = GermlineAmpDel.generateFilename(fileSources.Purple, sampleId);

            if(!Files.exists(Paths.get(germlineAmpDelFile)))
            {
                // try pre v3.0 germline deletions
                germlineAmpDelFile = germlineAmpDelFile.replaceAll("purple.germline_amp_del.tsv", "purple.germline.deletion.tsv");
            }

            List<GermlineAmpDel> germlineAmpDels = GermlineAmpDel.read(germlineAmpDelFile);
            germlineAmpDels.forEach(x -> comparableItems.add(new GermlineAmpDelData(x, mFields)));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to read germline amp-del data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }
}
