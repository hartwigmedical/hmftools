package com.hartwig.hmftools.lilac.qc;

import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;

public final class LilacQC
{
    private final Set<LilacQCStatus> mStatus;
    private final AminoAcidQC mAminoAcidQC;
    private final BamQC mBamQC;
    private final CoverageQC mCoverageQC;
    private final HaplotypeQC mHaplotypeQC;
    private final SomaticVariantQC mSomaticVariantQC;

    public final List<String> header()
    {
        return Lists.newArrayList();
        
//        return CollectionsKt.plus((Collection) CollectionsKt.plus((Collection) CollectionsKt.plus((Collection) CollectionsKt.plus((Collection) CollectionsKt
//                .plus((Collection) CollectionsKt.listOf((Object) "status"), (Iterable) aminoAcidQC.header()), (Iterable) bamQC.header()), (Iterable) coverageQC
//                .header()), (Iterable) haplotypeQC.header()), (Iterable) somaticVariantQC.header());
    }

    public final List<String> body()
    {
        return Lists.newArrayList();
//        return CollectionsKt.plus((Collection) CollectionsKt.plus((Collection) CollectionsKt.plus((Collection) CollectionsKt.plus((Collection) CollectionsKt
//                .plus((Collection) CollectionsKt.listOf((Object) CollectionsKt.joinToString$default((Iterable) status, (CharSequence) ",", null, null, (int) 0, null, null, (int) 62, null)), (Iterable) aminoAcidQC
//                        .body()), (Iterable) bamQC.body()), (Iterable) coverageQC.body()), (Iterable) haplotypeQC.body()), (Iterable) somaticVariantQC
//                .body());
    }

    public final void writefile(final String fileName)
    {
        /*
        Intrinsics.checkParameterIsNotNull((Object) fileName, (String) "fileName");
        String header =
                CollectionsKt.joinToString$default((Iterable) header(), (CharSequence) "\t", null, null, (int) 0, null, null, (int) 62, null);
        String body =
                CollectionsKt.joinToString$default((Iterable) body(), (CharSequence) "\t", null, null, (int) 0, null, null, (int) 62, null);
        File file = new File(fileName);
        FilesKt.writeText$default((File) file, (String) (header + "\n"), null, (int) 2, null);
        FilesKt.appendText$default((File) file, (String) (body + "\n"), null, (int) 2, null);
        
         */
    }


    public LilacQC(final Set<LilacQCStatus> status, final AminoAcidQC aminoAcidQC, final BamQC bamQC,
            final CoverageQC coverageQC, final HaplotypeQC haplotypeQC, final SomaticVariantQC somaticVariantQC)
    {
        mStatus = status;
        mAminoAcidQC = aminoAcidQC;
        mBamQC = bamQC;
        mCoverageQC = coverageQC;
        mHaplotypeQC = haplotypeQC;
        mSomaticVariantQC = somaticVariantQC;
    }


    public final LilacQC create(final AminoAcidQC aminoAcidQC, final BamQC bamQC, final CoverageQC coverageQC,
            final HaplotypeQC haplotypeQC, final SomaticVariantQC somaticVariantQC)
    {
        Set status = new LinkedHashSet();
        if(haplotypeQC.getUnusedHaplotypes() > 0)
        {
            status.add(LilacQCStatus.WARN_UNMATCHED_HAPLOTYPE);
        }
        if(coverageQC.ATypes == 0 || coverageQC.BTypes == 0 || coverageQC.CTypes == 0)
        {
            status.add(LilacQCStatus.WARN_UNMATCHED_TYPE);
        }
        if(bamQC.getDiscardedIndelFragments() > 0)
        {
            status.add(LilacQCStatus.WARN_UNMATCHED_INDEL);
        }
        if(somaticVariantQC.unmatchedVariants())
        {
            status.add(LilacQCStatus.WARN_UNMATCHED_SOMATIC_VARIANT);
        }
        if(coverageQC.getPercentWildcard() > 0.0)
        {
            status.add(LilacQCStatus.WARN_WILDCARD_MATCH);
        }
        if(status.isEmpty())
        {
            status.add(LilacQCStatus.PASS);
        }
        return new LilacQC(status, aminoAcidQC, bamQC, coverageQC, haplotypeQC, somaticVariantQC);
    }
}
