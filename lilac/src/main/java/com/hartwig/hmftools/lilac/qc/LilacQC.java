package com.hartwig.hmftools.lilac.qc;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.DELIM;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

public final class LilacQC
{
    private final Set<LilacQCStatus> mStatus;
    private final AminoAcidQC mAminoAcidQC;
    private final BamQC mBamQC;
    private final CoverageQC mCoverageQC;
    private final HaplotypeQC mHaplotypeQC;
    private final SomaticVariantQC mSomaticVariantQC;

    public final String header()
    {
        StringJoiner sj = new StringJoiner(DELIM);
        sj.add("status");
        mAminoAcidQC.header().forEach(x -> sj.add(x));
        mBamQC.header().forEach(x -> sj.add(x));
        mCoverageQC.header().forEach(x -> sj.add(x));
        mHaplotypeQC.header().forEach(x -> sj.add(x));
        mSomaticVariantQC.header().forEach(x -> sj.add(x));

        return sj.toString();
    }

    public final String body()
    {
        StringJoiner status = new StringJoiner(",");
        mStatus.forEach(x -> status.add(x.toString()));

        StringJoiner sj = new StringJoiner(DELIM);
        sj.add(status.toString());
        mAminoAcidQC.body().forEach(x -> sj.add(x));
        mBamQC.body().forEach(x -> sj.add(x));
        mCoverageQC.body().forEach(x -> sj.add(x));
        mHaplotypeQC.body().forEach(x -> sj.add(x));
        mSomaticVariantQC.body().forEach(x -> sj.add(x));

        return sj.toString();
    }

    public final void writefile(final String fileName)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write(header());
            writer.newLine();

            writer.write(body());
            writer.newLine();

            writer.close();
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to write {}: {}", fileName, e.toString());
            return;
        }
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


    public static LilacQC create(
            final AminoAcidQC aminoAcidQC, final BamQC bamQC, final CoverageQC coverageQC,
            final HaplotypeQC haplotypeQC, final SomaticVariantQC somaticVariantQC)
    {
        Set<LilacQCStatus> statusList = Sets.newHashSet();

        if(haplotypeQC.UnusedHaplotypes > 0)
        {
            statusList.add(LilacQCStatus.WARN_UNMATCHED_HAPLOTYPE);
        }

        if(coverageQC.ATypes == 0 || coverageQC.BTypes == 0 || coverageQC.CTypes == 0)
        {
            statusList.add(LilacQCStatus.WARN_UNMATCHED_TYPE);
        }

        if(bamQC.getDiscardedIndelFragments() > 0)
        {
            statusList.add(LilacQCStatus.WARN_UNMATCHED_INDEL);
        }

        if(somaticVariantQC.unmatchedVariants())
        {
            statusList.add(LilacQCStatus.WARN_UNMATCHED_SOMATIC_VARIANT);
        }

        if(coverageQC.PercentWildcard > 0.0)
        {
            statusList.add(LilacQCStatus.WARN_WILDCARD_MATCH);
        }

        if(statusList.isEmpty())
        {
            statusList.add(LilacQCStatus.PASS);
        }

        return new LilacQC(statusList, aminoAcidQC, bamQC, coverageQC, haplotypeQC, somaticVariantQC);
    }
}
