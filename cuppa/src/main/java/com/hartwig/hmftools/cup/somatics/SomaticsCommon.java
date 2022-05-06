package com.hartwig.hmftools.cup.somatics;

import static java.lang.Math.pow;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CupCalcs.convertToPercentages;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadRefSampleCounts;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Matrix;

public final class SomaticsCommon
{
    public static final String SPLIT_AID_APOBEC = "split_aid_apobec_gen_pos";
    public static final String SPLIT_AID_APOBEC_DESC = "Exclude 8 AID/APOBEC trinucleotide contexts from genomic positions";

    public static final String NORMALISE_COPY_NUMBER = "normalise_cn";
    public static final String NORMALISE_COPY_NUMBER_DESC = "Adjust genomic-position counts by copy number ";

    public static final String EXCLUDE_SNV_96_AID_APOBEC = "exclude_aid_apobec_snv_96";
    public static final String EXCLUDE_SNV_96_AID_APOBEC_DESC = "Exclude 8 AID/APOBEC trinucleotide contexts from SNV-96 counts";

    public static final String INCLUDE_AID_APOBEC_SIG = "aid_apobec_sig_feature";
    public static final String INCLUDE_AID_APOBEC_SIG_DESC = "Add an enriched AID/APOBEC signature feature";

    public static void applyMaxCssAdjustment(double maxCssScore, final Map<String,Double> cancerCssTotals, double adjustFactor)
    {
        if(adjustFactor == 0)
            return;

        double adjustedFactor = pow(maxCssScore, adjustFactor);

        for(Map.Entry<String,Double> entry : cancerCssTotals.entrySet())
        {
            double adjCancerScore = pow(entry.getValue(), adjustedFactor);
            cancerCssTotals.put(entry.getKey(), adjCancerScore);
        }

        convertToPercentages(cancerCssTotals);
    }

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
                return null;

            CUP_LOGGER.info("combined {} counts from {} samples", type, samplesList.size());

            final double[][] subData = subMatrix.getData();

            if(combinedMatrix == null)
            {
                combinedMatrix = new Matrix(subMatrix.Rows, refSampleCount);
            }

            for(int s = 0; s < samplesList.size(); ++s)
            {
                final String sampleId = samplesList.get(s);

                if(!refSampleIds.contains(sampleId))
                    continue;

                sampleCountsIndex.put(sampleId, sampleIndex);

                for(int r = 0; r < combinedMatrix.Rows; ++r)
                {
                    combinedMatrix.set(r, sampleIndex, subData[r][s]);
                }

                ++sampleIndex;
            }
        }

        combinedMatrix.cacheTranspose();

        return combinedMatrix;
    }

}
