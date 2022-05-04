package com.hartwig.hmftools.common.sigs;

import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;

import java.util.Map;

import com.hartwig.hmftools.common.variant.SomaticVariant;

public class SnvSigUtils
{
    public static final int SNV_TRINUCLEOTIDE_BUCKET_COUNT = 96;

    public static void populateBucketMap(final Map<String,Integer> bucketNameIndexMap)
    {
        char[] refBases = {'C', 'T'};
        char[] bases = {'A','C', 'G', 'T'};
        int index = 0;

        for(int i = 0; i < refBases.length; ++i)
        {
            char ref = refBases[i];

            for(int j = 0; j < bases.length; ++j)
            {
                char alt = bases[j];

                if(ref != alt)
                {
                    String baseChange = String.format("%c>%c", ref, alt);

                    for (int k = 0; k < bases.length; ++k)
                    {
                        char before = bases[k];

                        for (int l = 0; l < bases.length; ++l)
                        {
                            char after = bases[l];

                            String context = String.format("%c%c%c", before, ref, after);

                            String bucketName = baseChange + "_" + context;

                            bucketNameIndexMap.put(bucketName, index);
                            ++index;
                        }
                    }
                }
            }
        }
    }

    public static String contextFromVariant(final SomaticVariant variant)
    {
        return variantContext(variant.ref(), variant.alt(), variant.trinucleotideContext());
    }

    public static String variantContext(final String ref, final String alt, final String trinucContext)
    {
        // convert base change to standard set and the context accordingly
        if(ref.charAt(0) == 'A' || ref.charAt(0) == 'G')
        {
            return String.format("%c>%c_%c%c%c",
                    swapDnaBase(ref.charAt(0)), swapDnaBase(alt.charAt(0)),
                    swapDnaBase(trinucContext.charAt(2)), swapDnaBase(trinucContext.charAt(1)), swapDnaBase(trinucContext.charAt(0)));
        }
        else
        {
            return String.format("%c>%c_%s",ref.charAt(0), alt.charAt(0), trinucContext);
        }
    }
}
