package com.hartwig.hmftools.lilac.qc;

import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.DELIM;

import com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage;
import com.hartwig.hmftools.lilac.coverage.HlaComplexCoverage;
import com.hartwig.hmftools.lilac.variant.SomaticCodingCount;
import com.hartwig.hmftools.lilac.hla.HlaCopyNumber;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public class HlaOut
{
    private final HlaComplexCoverage mReferenceCoverage;
    private final HlaComplexCoverage mTumorCoverage;
    private final List<HlaCopyNumber> mTumorCopyNumber;
    private final List<SomaticCodingCount> mSomaticCodingCount;

    public HlaOut(
            final HlaComplexCoverage referenceCoverage, final HlaComplexCoverage tumorCoverage,
            final List<HlaCopyNumber> tumorCopyNumber, final List<SomaticCodingCount> somaticCodingCount)
    {
        mReferenceCoverage = referenceCoverage;
        mTumorCoverage = tumorCoverage;
        mTumorCopyNumber = tumorCopyNumber;
        mSomaticCodingCount = somaticCodingCount;
    }

    public final void write(final String fileName)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write(generateAlleleHeader());
            writer.newLine();

            for(int i = 0; i < 6; ++i)
            {
                writer.write(generateAlleleBody(i));
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to write {}: {}", fileName, e.toString());
            return;
        }
    }

    private final String generateAlleleHeader()
    {
        StringJoiner header = new StringJoiner(DELIM)
                .add("Allele")
                .add("RefTotal")
                .add("RefUnique")
                .add("RefShared")
                .add("RefWild")
                .add("TumorTotal")
                .add("TumorUnique")
                .add("TumorShared")
                .add("TumorWild")
                .add("TumorCopyNumber")
                .add("SomaticMissense")
                .add("SomaticNonsenseOrFrameshift")
                .add("SomaticSplice")
                .add("SomaticSynonymous")
                .add("SomaticInframeIndel");
        return header.toString();
    }

    private final String generateAlleleBody(int index)
    {
        HlaAlleleCoverage ref = mReferenceCoverage.getAlleleCoverage().get(index);

        HlaAlleleCoverage tumor = !mTumorCoverage.getAlleleCoverage().isEmpty() ?
                mTumorCoverage.getAlleleCoverage().get(index) : new HlaAlleleCoverage(ref.Allele, 0, 0, 0);

        double copyNumber = mTumorCopyNumber.get(index).CopyNumber;
        SomaticCodingCount codingCount = mSomaticCodingCount.get(index);

        StringJoiner header = new StringJoiner(DELIM)
                .add(ref.Allele.asFourDigit().toString())
                .add(String.valueOf(round(ref.TotalCoverage)))
                .add(String.valueOf(ref.UniqueCoverage))
                .add(String.valueOf(round(ref.SharedCoverage)))
                .add(String.valueOf(round(ref.WildCoverage)))
                .add(String.valueOf(round(tumor.TotalCoverage)))
                .add(String.valueOf(tumor.UniqueCoverage))
                .add(String.valueOf(round(tumor.SharedCoverage)))
                .add(String.valueOf(round(tumor.WildCoverage)))
                .add(String.format("%.2f", copyNumber))
                .add(String.format("%.2f", codingCount.Missense))
                .add(String.format("%.2f", codingCount.Nonsense))
                .add(String.format("%.2f", codingCount.Splice))
                .add(String.format("%.2f", codingCount.Synonymous))
                .add(String.format("%.2f", codingCount.InframeIndel));

        return header.toString();
    }

    public static HlaOut create(
            final HlaComplexCoverage referenceCoverage, final HlaComplexCoverage tumorCoverage,
            final List<HlaCopyNumber> tumorCopyNumber, final List<SomaticCodingCount> somaticCodingCount)
    {
        List<HlaCopyNumber> sortedCopyNumber = tumorCopyNumber.stream().collect(Collectors.toList());
        Collections.sort(sortedCopyNumber, new HlaCopyNumberSorter());

        List<SomaticCodingCount> sortedCodingCount = somaticCodingCount.stream().collect(Collectors.toList());
        Collections.sort(sortedCodingCount, new SomaticCodingCountSorter());

        return new HlaOut(referenceCoverage, tumorCoverage, sortedCopyNumber, sortedCodingCount);
    }

    private static class HlaCopyNumberSorter implements Comparator<HlaCopyNumber>
    {
        public int compare(final HlaCopyNumber first, final HlaCopyNumber second)
        {
            return first.Allele.compareTo(second.Allele);
        }
    }

    private static class SomaticCodingCountSorter implements Comparator<SomaticCodingCount>
    {
        public int compare(final SomaticCodingCount first, final SomaticCodingCount second)
        {
            return first.Allele.compareTo(second.Allele);
        }
    }

}
