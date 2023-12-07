package com.hartwig.hmftools.common.cuppa2;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

public class Cuppa2DataFiles
{
    public static String generatePredictionsPath(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + ".cuppa_vis_data.tsv";
    }

    public static String generateVisualizationPath(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + ".cuppa_vis.tsv";
    }

    public static String generatePredSummPath(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + ".pred_summ.tsv";
    }
}
