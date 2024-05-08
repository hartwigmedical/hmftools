package com.hartwig.hmftools.common.variant;

import java.util.List;

import htsjdk.variant.vcf.VCFHeader;

public class GenotypeIds
{
    public final int ReferenceOrdinal;
    public final int TumorOrdinal;
    public final String ReferenceId;
    public final String TumorId;

    public GenotypeIds(final int referenceOrdinal, final int tumorOrdinal, final String referenceId, final String tumorId)
    {
        ReferenceOrdinal = referenceOrdinal;
        TumorOrdinal = tumorOrdinal;
        ReferenceId = referenceId;
        TumorId = tumorId;
    }

    public boolean hasReference() { return ReferenceOrdinal >= 0; }
    public boolean hasTumor() { return TumorOrdinal >= 0; }

    public static boolean hasValidSampleIds(
            final VCFHeader header, final String referenceId, final String tumorId, boolean referenceFirst, boolean allowContains)
    {
        // ordering is only relevant if both are supplied
        List<String> vcfSampleNames = header.getGenotypeSamples();

        int expectedRefOrdinal;
        int expectedTumorOrdinal;

        if(tumorId != null && !tumorId.isEmpty() && referenceId != null && !referenceId.isEmpty())
        {
            expectedRefOrdinal = referenceFirst ? 0 : 1;
            expectedTumorOrdinal = referenceFirst ? 1 : 0;
        }
        else if(tumorId != null && !tumorId.isEmpty())
        {
            expectedRefOrdinal = -1;
            expectedTumorOrdinal = 0;
        }
        else
        {
            expectedRefOrdinal = 0;
            expectedTumorOrdinal = -1;
        }

        for(int i = 0; i < vcfSampleNames.size(); ++i)
        {
            String vcfSampleName = vcfSampleNames.get(i);

            if(vcfSampleName.equals(referenceId) || (allowContains && !referenceId.isEmpty() && vcfSampleName.contains(referenceId)))
            {
                if(expectedRefOrdinal >= 0 && i != expectedRefOrdinal)
                    return false;
            }
            else if(vcfSampleName.equals(tumorId) || (allowContains && !tumorId.isEmpty() && vcfSampleName.contains(tumorId)))
            {
                if(expectedTumorOrdinal >= 0 && i != expectedTumorOrdinal)
                    return false;
            }
        }

        return true;
    }

    public static GenotypeIds fromVcfHeader(final VCFHeader header, final String referenceId, final String tumorId)
    {
        List<String> vcfSampleNames = header.getGenotypeSamples();

        int tumorOrdinal = -1;
        int referenceOrdinal = -1;
        String vcfTumorId = "";
        String vcfRefefenceId = "";

        for(int test = 0; test <= 1; ++test)
        {
            boolean checkExact = test == 0;

            for(int i = 0; i < vcfSampleNames.size(); ++i)
            {
                String vcfSampleName = vcfSampleNames.get(i);

                if(tumorOrdinal < 0 && tumorId != null && (vcfSampleName.equals(tumorId)
                || (!checkExact && tumorId != null && vcfSampleName.contains(tumorId))))
                {
                    vcfTumorId = vcfSampleNames.get(i);
                    tumorOrdinal = i;
                }
                else if(referenceOrdinal < 0 && referenceId != null && (vcfSampleName.equals(referenceId)
                || (!checkExact && referenceId != null && vcfSampleName.contains(referenceId))))
                {
                    vcfRefefenceId = vcfSampleNames.get(i);
                    referenceOrdinal = i;
                }
            }
        }

        return new GenotypeIds(referenceOrdinal, tumorOrdinal, vcfRefefenceId, vcfTumorId);
    }
}
