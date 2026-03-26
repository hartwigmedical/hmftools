package com.hartwig.hmftools.common.pipeline;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.enforceChrPrefix;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

public final class MiscToolFiles
{
    public static String SAGE_VIS_PLOT_FILE_EXTENSION = ".html.gz";

    public static String generateSageVisFilePrefix(
            final String geneName, final String chromosome, final int position, final String ref, final String alt)
    {
        // of the form .sage.GENE_chr1_1000_A_AC
        return format(".sage.%s_%s_%d_%s_%s", geneName, enforceChrPrefix(chromosome), position, ref, alt);
    }

    public static String QSEE_PLOT_FILE_ID = ".qsee.vis.report.png";

    public static String generateQseeVisPlot(final String basePath, final String sampleId)
    {
        return checkAddDirSeparator(basePath) + sampleId + QSEE_PLOT_FILE_ID;
    }

}
