package com.hartwig.hmftools.cup.somatics;

import static java.lang.Math.pow;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CupCalcs.convertToPercentages;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadRefSampleCounts;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Matrix;

public final class SomaticsCommon
{
    public static final String INCLUDE_AID_APOBEC = "include_aid_apobec_gen_pos";
    public static final String INCLUDE_AID_APOBEC_DESC = "Include 8 AID/APOBEC trinucleotide contexts in genomic positions";

    public static final String INCLUDE_AID_APOBEC_SIG = "aid_apobec_sig_feature";
    public static final String INCLUDE_AID_APOBEC_SIG_DESC = "Add an enriched AID/APOBEC signature feature";

    public static final String EXCLUDE_GEN_POS_CHR_X = "exclude_gen_pos_chr_x";
    public static final String EXCLUDE_GEN_POS_CHR_X_DESC = "Exclude chromosome X from genomic position";

    public static final String GEN_POS_MAX_SAMPLE_COUNT_CFG = "gen_pos_max_sample_count";
    public static final String GEN_POS_MAX_SAMPLE_COUNT_DESC = "Max SNV count per sample for genomic position";

    public static final String GEN_POS_BUCKET_SIZE_CFG = "gen_pos_bucket_size";
    public static final String GEN_POS_BUCKET_SIZE_DESC = "Genomic position bucket size (default: 500K)";

    public static final String INTEGER_FORMAT = "%.0f";
    public static final String DEC_1_FORMAT = "%.1f";
    public static final String DEC_3_FORMAT = "%.3f";

    public static Matrix loadMultipleMatrixFiles(
            final List<String> filenames, final List<String> refSampleIds, final Map<String,Integer> sampleCountsIndex, final String type)
    {
        int refSampleCount = refSampleIds.size();

        Matrix combinedMatrix = null;
        int sampleIndex = 0;

        for(String filename : filenames)
        {
            final List<String> samplesList = Lists.newArrayList();
            final Matrix subMatrix = loadRefSampleCounts(filename, samplesList, Lists.newArrayList("BucketName"));

            if(subMatrix == null)
            {
                CUP_LOGGER.info("failed to load type({}) file{{}}", type, filename);
                return null;
            }

            CUP_LOGGER.info("combined {} counts from {} samples", type, samplesList.size());

            final double[][] subData = subMatrix.getData();

            if(combinedMatrix == null)
            {
                combinedMatrix = new Matrix(refSampleCount, subMatrix.Cols);
            }

            final double[][] combinedData = combinedMatrix.getData();

            for(int s = 0; s < samplesList.size(); ++s)
            {
                final String sampleId = samplesList.get(s);

                if(!refSampleIds.contains(sampleId))
                    continue;

                sampleCountsIndex.put(sampleId, sampleIndex);

                for(int b = 0; b < combinedMatrix.Cols; ++b)
                {
                    combinedData[sampleIndex][b] = subData[s][b];
                }

                ++sampleIndex;
            }
        }

        return combinedMatrix;
    }
}
