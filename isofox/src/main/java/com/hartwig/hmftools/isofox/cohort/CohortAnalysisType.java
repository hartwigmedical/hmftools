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
    GENE_DISTRIBUTION, // creates gene expression percentiles by cancer type and pan-cancer
    TRANSCRIPT_DISTRIBUTION, // as above but for transcripts
    SAMPLE_GENE_PERCENTILES, // reports where on a percentile a specific sample's gene expression lies
    ALT_SPLICE_JUNCTION, // combine and analyse alternate splice junctions for a cohort
    SPLICE_VARIANT_MATCHING, // match alternate splice junctions with (candidate-splicing) somatic variants
    FUSION, // process fusions for a cohort - filter passing fusions, form a cohort file, compare with external fusions
    PASSING_FUSION,
    RETAINED_INTRON,
    SAMPLE_ROUTINES,
    EXTERNAL_EXPRESSION_COMPARE, // combine expression data from Isofox and another source
    GENE_EXPRESSION_COMPARE, // compare gene expression across 2 cohorts of samples
    GENE_EXPRESSION_MATRIX, // generates a single gene by sample matrix for expression data
    TRANSCRIPT_EXPRESSION_MATRIX; // as above but for transcript expression

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
