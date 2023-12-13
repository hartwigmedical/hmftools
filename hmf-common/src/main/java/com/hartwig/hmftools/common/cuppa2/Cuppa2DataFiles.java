package com.hartwig.hmftools.common.cuppa2;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

public class Cuppa2DataFiles
{
    public static String generateVisDataPath(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + ".cuppa.vis_data.tsv";
    }

    public static String generateVisPlotPath(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + ".cuppa.vis.png";
    }

    public static String generatePredSummPath(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + ".cuppa.pred_summ.tsv";
    }
}
