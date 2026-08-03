package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;

import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public abstract class ImageComparer extends ItemComparer
{
    public static final String FLD_DIMENSIONS = "Dimensions";
    public static final String FLD_PIXELS = "Pixels";

    private final Double mPixelAbsoluteThreshold;
    private final Double mPixelPercentThreshold;

    protected ImageComparer(
            final ComparConfig config, final Double pixelAbsoluteThreshold, final Double pixelPercentThreshold,
            final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mPixelAbsoluteThreshold = pixelAbsoluteThreshold;
        mPixelPercentThreshold = pixelPercentThreshold;

        mFields.add(new FieldInfo(FLD_DIMENSIONS, getOrMakeFieldCheck(fieldCheckMap, FLD_DIMENSIONS), null));

        mFields.add(new FieldInfo(
                FLD_PIXELS,
                getOrMakeFieldCheck(fieldCheckMap, FLD_PIXELS, mPixelAbsoluteThreshold, mPixelPercentThreshold),
                "%.2f"));
    }

    /*
    @VisibleForTesting
    public static List<Field> buildFields(final Double pixelAbsoluteThreshold, final Double pixelPercentThreshold)
    {
        return List.of(
                new StringField(FLD_DIMENSIONS, i -> ((ComparableImage) i).dimensionString(), true),
                new PixelField(FLD_PIXELS, i -> ((ComparableImage) i).Image, true,
                        pixelAbsoluteThreshold, pixelPercentThreshold)
        );
    }
    */
}
