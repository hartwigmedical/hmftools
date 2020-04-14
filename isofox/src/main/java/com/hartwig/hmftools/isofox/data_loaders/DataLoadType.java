package com.hartwig.hmftools.isofox.data_loaders;

import static com.hartwig.hmftools.isofox.novel.AltSpliceJunctionFinder.ALT_SJ_FILE_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.GENE_RESULTS_FILE;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.SUMMARY_FILE;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.TRANSCRIPT_RESULTS_FILE;

public enum DataLoadType
{
    SUMMARY,
    GENE_DISTRIBUTION,
    SAMPLE_GENE_PERCENTILES,
    TRANSCRIPT_DISTRIBUTION,
    ALT_SPLICE_JUNCTION,
    RETAINED_INTRON;

    public static String getFileId(DataLoadType type)
    {
        switch(type)
        {
            case GENE_DISTRIBUTION:
            case SAMPLE_GENE_PERCENTILES:
                return GENE_RESULTS_FILE;

            case TRANSCRIPT_DISTRIBUTION:
                return TRANSCRIPT_RESULTS_FILE;

            case SUMMARY:
                return SUMMARY_FILE;

            case ALT_SPLICE_JUNCTION:
                return ALT_SJ_FILE_ID;

            default:
                return "";
        }
    }
}
