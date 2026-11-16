package com.hartwig.hmftools.compar.cuppa;

import static com.hartwig.hmftools.compar.common.CategoryType.CUPPA_IMAGE;

import java.util.Collections;

import com.hartwig.hmftools.compar.common.ComparableImage;
import com.hartwig.hmftools.compar.common.CategoryType;

public class CuppaImageData extends ComparableImage
{
    public CuppaImageData(final String name, final String path)
    {
        super(name, path, Collections.emptyList());

        // PixelFieldValue
    }

    @Override
    public CategoryType category() { return CUPPA_IMAGE; }
}
