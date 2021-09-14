package com.hartwig.hmftools.isofox.cohort;

import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.ALT_SJ_FILE_ID;
import static com.hartwig.hmftools.common.rna.GeneExpressionFile.GENE_DATA_FILE_ID;
import static com.hartwig.hmftools.isofox.fusion.FusionWriter.RAW_FUSION_FILE_ID;
import static com.hartwig.hmftools.isofox.fusion.FusionWriter.PASS_FUSION_FILE_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.SPLICE_SITE_FILE;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.SUMMARY_FILE;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.TRANSCRIPT_RESULTS_FILE;

public enum AnalysisType
{
    SUMMARY,
    EXPRESSION_DISTRIBUTION, // produce pan-cancer and per-cancer median and percentile expression data
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
            case GENE_EXPRESSION_COMPARE:
                return GENE_DATA_FILE_ID;

            case SUMMARY:
                return SUMMARY_FILE;

            case ALT_SPLICE_JUNCTION:
            case SPLICE_VARIANT_MATCHING:
            case ALT_SPLICE_JUNCTION_MATRIX:
                return ALT_SJ_FILE_ID;

            case SPLICE_SITE_PERCENTILES:
                return SPLICE_SITE_FILE;

            case FUSION:
                return RAW_FUSION_FILE_ID;

            case PASSING_FUSION:
                return PASS_FUSION_FILE_ID;

            case GENE_EXPRESSION_MATRIX:
                return GENE_DATA_FILE_ID;

            case TRANSCRIPT_EXPRESSION_MATRIX:
                return TRANSCRIPT_RESULTS_FILE;

            default:
                return "";
        }
    }
}
