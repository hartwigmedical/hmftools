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
        // convert base change to standard set and the context accordingly
        String baseChange;
        String context;
        final String rawContext = variant.trinucleotideContext();

        if(variant.ref().charAt(0) == 'A' || variant.ref().charAt(0) == 'G')
        {
            baseChange = String.format("%c>%c", swapDnaBase(variant.ref().charAt(0)), swapDnaBase(variant.alt().charAt(0)));

            // convert the context as well
            context = String.format("%c%c%c",
                    swapDnaBase(rawContext.charAt(2)), swapDnaBase(rawContext.charAt(1)), swapDnaBase(rawContext.charAt(0)));
        }
        else
        {
            baseChange = variant.ref() + ">" + variant.alt();
            context = rawContext;
        }

        return baseChange + "_" + context;
    }
}
