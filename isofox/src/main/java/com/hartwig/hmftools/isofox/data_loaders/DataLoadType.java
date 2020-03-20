package com.hartwig.hmftools.isofox.data_loaders;

import static com.hartwig.hmftools.isofox.results.ResultsWriter.GENE_RESULTS_FILE;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.SUMMARY_FILE;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.TRANSCRIPT_RESULTS_FILE;

public enum DataLoadType
{
    SUMMARY,
    GENE,
    TRANSCRIPT,
    ALT_SPLICE_JUNCTIONS,
    RETAINED_INTRONS;

    public static String getFileId(DataLoadType type)
    {
        switch(type)
        {
            case GENE: return GENE_RESULTS_FILE;
            case TRANSCRIPT: return TRANSCRIPT_RESULTS_FILE;
            case SUMMARY: return SUMMARY_FILE;
            default: return "";
        }
    }
}
