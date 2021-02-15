package com.hartwig.hmftools.isofox.cohort;

import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.ALT_SJ_FILE_ID;
import static com.hartwig.hmftools.isofox.fusion.FusionWriter.FUSION_FILE_ID;
import static com.hartwig.hmftools.isofox.fusion.cohort.FusionCohort.PASS_FUSION_FILE_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.GENE_RESULTS_FILE;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.SPLICE_SITE_FILE;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.SUMMARY_FILE;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.TRANSCRIPT_RESULTS_FILE;

public enum AnalysisType
{
    SUMMARY,
    GENE_DISTRIBUTION, // creates gene expression percentiles by cancer type and pan-cancer
    TRANSCRIPT_DISTRIBUTION, // as above but for transcripts
    EXPRESSION_COHORT_MEDIANS, // produce pan-cancer and per-cancer median expression values
    SAMPLE_GENE_PERCENTILES, // reports where on a percentile a specific sample's gene expression lies
    ALT_SPLICE_JUNCTION, // combine and analyse alternate splice junctions for a cohort
    ALT_SPLICE_JUNCTION_MATRIX, // generates a for a cohort's alt-SJs
    SPLICE_VARIANT_MATCHING, // match alternate splice junctions with (candidate-splicing) somatic variants
    RECURRENT_SPLICE_VARIANTS, // find recurrent splice variants, input is a somatic table file
    FUSION, // process fusions for a cohort - filter passing fusions, form a cohort file, compare with external fusions
    PASSING_FUSION, // produce per-sample filtered fusions
    RETAINED_INTRON, // produce cohort file for retained introns
    SPLICE_SITE_PERCENTILES, // produce cohort percentage-spiced-in data from per-sample splice data files
    EXTERNAL_EXPRESSION_COMPARE, // combine expression data from Isofox and another source
    GENE_EXPRESSION_COMPARE, // compare gene expression across 2 cohorts of samples
    GENE_EXPRESSION_MATRIX, // generates a matrix for gene expression data
    TRANSCRIPT_EXPRESSION_MATRIX; // as above but for transcript expression

    public static String getIsofoxFileId(AnalysisType type)
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
            case ALT_SPLICE_JUNCTION_MATRIX:
                return ALT_SJ_FILE_ID;

            case SPLICE_SITE_PERCENTILES:
                return SPLICE_SITE_FILE;

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
