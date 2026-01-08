package com.hartwig.hmftools.compar.cuppa;

import static com.hartwig.hmftools.compar.common.Category.CUPPA_IMAGE;

import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.ComparableImage;

public class CuppaImageData extends ComparableImage
{
    public static final String FLD_VIS_IMAGE = "cuppa_vis_image";

    public CuppaImageData(final String name, final String path)
    {
        super(name, path);
    }

    @Override
    public Category category() { return CUPPA_IMAGE; }
}
