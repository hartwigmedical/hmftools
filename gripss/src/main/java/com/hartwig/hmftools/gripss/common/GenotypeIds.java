package com.hartwig.hmftools.gripss.common;

import static com.hartwig.hmftools.gripss.GripssConfig.GR_LOGGER;

import java.util.List;

import htsjdk.variant.vcf.VCFHeader;

public class GenotypeIds
{
    public final int ReferenceOrdinal;
    public final int TumorOrdinal;
    public final String ReferenceId;
    public final String TumorId;
    public final boolean GermlineMode;

    public GenotypeIds(final int referenceOrdinal, final int tumorOrdinal, final String referenceId, final String tumorId, final boolean germlineMode)
    {
        ReferenceOrdinal = referenceOrdinal;
        TumorOrdinal = tumorOrdinal;
        ReferenceId = referenceId;
        TumorId = tumorId;
        GermlineMode = germlineMode;
    }

    public boolean hasReference() { return ReferenceOrdinal >= 0; }

    public static GenotypeIds fromVcfHeader(
            final VCFHeader header, final String referenceId, final String tumorId, final boolean germlineMode)
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

                if(tumorOrdinal < 0 && (vcfSampleName.equals(tumorId) || (!checkExact && vcfSampleName.contains(tumorId))))
                {
                    vcfTumorId = vcfSampleNames.get(i);
                    tumorOrdinal = i;
                }
                else if(referenceOrdinal < 0 && (vcfSampleName.equals(referenceId) || (!checkExact && vcfSampleName.contains(referenceId))))
                {
                    vcfRefefenceId = vcfSampleNames.get(i);
                    referenceOrdinal = i;
                }
            }
        }

        if(tumorOrdinal < 0 || (!referenceId.isEmpty() && referenceOrdinal < 0))
        {
            GR_LOGGER.error("missing sample names in VCF: {}", vcfSampleNames);
            return null;
        }

        return new GenotypeIds(referenceOrdinal, tumorOrdinal, vcfRefefenceId, vcfTumorId, germlineMode);
    }
}
