package com.hartwig.hmftools.compar;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.field.Field;
import com.hartwig.hmftools.compar.common.field.ImageField;
import com.hartwig.hmftools.compar.common.field.StringField;

public abstract class ImageComparer implements ItemComparer
{
    public static final String FLD_DIMENSIONS = "Dimensions";
    public static final String FLD_PIXELS = "Pixels";

    private final Double mPixelAbsoluteThreshold;
    private final Double mPixelPercentThreshold;

    protected ImageComparer(final Double pixelAbsoluteThreshold, final Double pixelPercentThreshold)
    {
        mPixelAbsoluteThreshold = pixelAbsoluteThreshold;
        mPixelPercentThreshold = pixelPercentThreshold;
    }

    @Override
    public List<Field> fields(final MatchLevel matchLevel)
    {
        return buildFields(mPixelAbsoluteThreshold, mPixelPercentThreshold);
    }

    @VisibleForTesting
    public static List<Field> buildFields(final Double pixelAbsoluteThreshold, final Double pixelPercentThreshold)
    {
        return List.of(
                new StringField(FLD_DIMENSIONS, i -> ((ComparableImage) i).dimensionString(), true),
                new ImageField(FLD_PIXELS, i -> ((ComparableImage) i).Image, true, pixelAbsoluteThreshold, pixelPercentThreshold)
        );
    }
}
