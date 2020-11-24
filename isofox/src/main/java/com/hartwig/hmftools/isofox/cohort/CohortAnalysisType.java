package com.hartwig.hmftools.isofox.cohort;

import static com.hartwig.hmftools.isofox.fusion.FusionWriter.FUSION_FILE_ID;
import static com.hartwig.hmftools.isofox.fusion.cohort.FusionCohort.PASS_FUSION_FILE_ID;
import static com.hartwig.hmftools.isofox.novel.AltSpliceJunctionFinder.ALT_SJ_FILE_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.GENE_RESULTS_FILE;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.SUMMARY_FILE;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.TRANSCRIPT_RESULTS_FILE;

public enum CohortAnalysisType
{
    SUMMARY,
    GENE_DISTRIBUTION,
    SAMPLE_GENE_PERCENTILES,
    TRANSCRIPT_DISTRIBUTION,
    ALT_SPLICE_JUNCTION,
    SPLICE_VARIANT_MATCHING,
    FUSION,
    PASSING_FUSION,
    RETAINED_INTRON,
    SAMPLE_ROUTINES,
    EXTERNAL_EXPRESSION_COMPARE,
    GENE_EXPRESSION_COMPARE,
    GENE_EXPRESSION_MATRIX,
    TRANSCRIPT_EXPRESSION_MATRIX;

    public static String getFileId(CohortAnalysisType type)
    {
        switch(type)
        {
            case GENE_DISTRIBUTION:
            case SAMPLE_GENE_PERCENTILES:
            case GENE_EXPRESSION_COMPARE:
                return GENE_RESULTS_FILE;

            case TRANSCRIPT_DISTRIBUTION:
                return TRANSCRIPT_RESULTS_FILE;

            case SUMMARY:
                return SUMMARY_FILE;

            case ALT_SPLICE_JUNCTION:
            case SPLICE_VARIANT_MATCHING:
                return ALT_SJ_FILE_ID;

            case FUSION:
                return FUSION_FILE_ID;

            case PASSING_FUSION:
                return PASS_FUSION_FILE_ID;

            case GENE_EXPRESSION_MATRIX:
                return GENE_RESULTS_FILE;

            case TRANSCRIPT_EXPRESSION_MATRIX:
                return TRANSCRIPT_RESULTS_FILE;

            default:
                return "";
        }
    }
}
