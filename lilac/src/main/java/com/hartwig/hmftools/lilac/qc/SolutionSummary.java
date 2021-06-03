package com.hartwig.hmftools.lilac.qc;

import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.DELIM;

import com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage;
import com.hartwig.hmftools.lilac.coverage.HlaComplexCoverage;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.variant.SomaticCodingCount;
import com.hartwig.hmftools.lilac.variant.CopyNumberAssignment;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public class SolutionSummary
{
    private final HlaComplexCoverage mReferenceCoverage;
    private final HlaComplexCoverage mTumorCoverage;
    private final Map<HlaAllele,Double> mTumorCopyNumber;
    private final List<SomaticCodingCount> mSomaticCodingCount;

    public SolutionSummary(
            final HlaComplexCoverage referenceCoverage, final HlaComplexCoverage tumorCoverage,
            final Map<HlaAllele,Double> tumorCopyNumber, final List<SomaticCodingCount> somaticCodingCount)
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

        double copyNumber = mTumorCopyNumber.get(ref.Allele);
        SomaticCodingCount codingCount = mSomaticCodingCount.get(index);

        StringJoiner header = new StringJoiner(DELIM)
                .add(ref.Allele.toString())
                .add(String.valueOf(round(ref.TotalCoverage)))
                .add(String.valueOf(ref.UniqueCoverage))
                .add(String.valueOf(round(ref.SharedCoverage)))
                .add(String.valueOf(round(ref.WildCoverage)))
                .add(String.valueOf(round(tumor.TotalCoverage)))
                .add(String.valueOf(tumor.UniqueCoverage))
                .add(String.valueOf(round(tumor.SharedCoverage)))
                .add(String.valueOf(round(tumor.WildCoverage)))
                .add(String.format("%.2f", copyNumber))
                .add(String.format("%.2f", codingCount.missense()))
                .add(String.format("%.2f", codingCount.nonsense()))
                .add(String.format("%.2f", codingCount.splice()))
                .add(String.format("%.2f", codingCount.synonymous()))
                .add(String.format("%.2f", codingCount.inframeIndel()));

        return header.toString();
    }

    public static SolutionSummary create(
            final HlaComplexCoverage referenceCoverage, final HlaComplexCoverage tumorCoverage,
            final Map<HlaAllele,Double> tumorCopyNumber, final List<SomaticCodingCount> somaticCodingCount)
    {
        List<SomaticCodingCount> sortedCodingCount = somaticCodingCount.stream().collect(Collectors.toList());
        Collections.sort(sortedCodingCount, new SomaticCodingCountSorter());

        return new SolutionSummary(referenceCoverage, tumorCoverage, tumorCopyNumber, sortedCodingCount);
    }

    private static class SomaticCodingCountSorter implements Comparator<SomaticCodingCount>
    {
        public int compare(final SomaticCodingCount first, final SomaticCodingCount second)
        {
            return first.Allele.compareTo(second.Allele);
        }
    }

}
