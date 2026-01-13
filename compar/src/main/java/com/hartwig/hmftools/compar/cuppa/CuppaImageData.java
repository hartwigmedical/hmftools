package com.hartwig.hmftools.compar.cuppa;

import static com.hartwig.hmftools.compar.common.CategoryType.CUPPA_IMAGE;

import com.hartwig.hmftools.compar.ComparableImage;
import com.hartwig.hmftools.compar.common.CategoryType;

public class CuppaImageData extends ComparableImage
{
    public CuppaImageData(final String name, final String path)
    {
        super(name, path);
    }

    @Override
    public CategoryType category() { return CUPPA_IMAGE; }
}
