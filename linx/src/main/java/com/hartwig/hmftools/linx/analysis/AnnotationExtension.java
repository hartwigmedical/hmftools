package com.hartwig.hmftools.linx.analysis;

import static com.hartwig.hmftools.linx.LinxOutput.ITEM_DELIM;

import java.util.List;

import com.google.common.collect.Lists;

public enum AnnotationExtension
{
    DOUBLE_MINUTES,
    CANDIDATE_VIS_DOUBLE_MINUTES,
    LINE_CHAINS,
    UNDER_CLUSTERING;

    public static List<AnnotationExtension> fromConfig(final String configStr)
    {
        List<AnnotationExtension> extensions = Lists.newArrayList();

        if(configStr == null || configStr.isEmpty())
            return extensions;

        for(String extension : configStr.split(ITEM_DELIM, -1))
        {
            if(extension.equals("DM")) // old value
                extensions.add(DOUBLE_MINUTES);
            else
                extensions.add(AnnotationExtension.valueOf(extension));
        }

        return extensions;
    }
}
