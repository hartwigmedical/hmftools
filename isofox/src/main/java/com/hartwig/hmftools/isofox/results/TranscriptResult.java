package com.hartwig.hmftools.isofox.results;

import static com.hartwig.hmftools.isofox.common.GeneCollection.TRANS_COUNT;
import static com.hartwig.hmftools.isofox.common.GeneCollection.UNIQUE_TRANS_COUNT;
import static com.hartwig.hmftools.isofox.common.RegionReadData.findExonRegion;
import static com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesGenerator.FL_FREQUENCY;
import static com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesGenerator.FL_LENGTH;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_NAME;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_TRANS_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_TRANS_NAME;

import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.isofox.common.FragmentMatchType;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.RegionReadData;

import org.immutables.value.Value;

@Value.Immutable
public abstract class TranscriptResult
{
    public abstract TranscriptData trans();

    public abstract int exonsFound();
    public abstract int shortSupportingFragments();
    public abstract int shortUniqueFragments();
    public abstract int longSupportingFragments();
    public abstract int longUniqueFragments();
    public abstract int spliceJunctionsSupported();
    public abstract int spliceJunctionFragments();
    public abstract int spliceJunctionUniqueFragments();
    public abstract int uniqueSpliceJunctions();
    public abstract int uniqueSpliceJunctionsSupported();
    public abstract int exonicBases();
    public abstract int exonicBaseCoverage();
    public abstract int uniqueBases();
    public abstract int uniqueBaseCoverage();
    public abstract double uniqueBaseAvgDepth();

    public abstract double effectiveLength();

    private double mFitAllocation;
    private double mTPM;

    public static TranscriptResult createTranscriptResults(
            final GeneCollection geneCollection, final GeneReadData geneReadData, final TranscriptData transData,
            final List<int[]> expRateFragmentLengths)
    {
        int exonsFound = 0;

        int spliceJunctionsSupported = 0;
        int exonicBases = 0;
        int exonicBaseCoverage = 0;

        int uniqueExonicBases = 0;
        int uniqueExonicBaseCoverage = 0;
        int uniqueExonicBaseTotalDepth = 0;

        int uniqueSpliceJunctions = 0;
        int uniqueSpliceJunctionsSupported = 0;

        /* Criteria for transcript selection
        - all exon junctions covered
        - unique exon junctions
        - split reads skipping exons
        - unique exon reads (but could cover introns as well
         */

        final List<ExonData> exons = transData.exons();

        for(int i = 0; i < exons.size(); ++i)
        {
            ExonData exon = exons.get(i);

            final RegionReadData exonReadData = findExonRegion(geneReadData.getExonRegions(), exon.ExonStart, exon.ExonEnd);
            if(exonReadData == null)
                continue;

            int exonCoverage = exonReadData.baseCoverage(1);
            exonicBaseCoverage += exonCoverage;

            if(exonCoverage > 0)
                ++exonsFound;

            if(exonReadData.getTransExonRefs().size() == 1)
            {
                uniqueExonicBases += exonReadData.uniqueBaseCount();
                uniqueExonicBaseCoverage += exonReadData.uniqueBaseCoverage(1);
                uniqueExonicBaseTotalDepth += exonReadData.uniqueBaseTotalDepth();
            }

            exonicBases += exon.ExonEnd - exon.ExonStart + 1;

            if(i > 0)
            {
                int[] sjReads = exonReadData.getTranscriptJunctionMatchCount(transData.TransId, SE_START);

                final ExonData prevExon = exons.get(i - 1);
                boolean sjUnique = isSpliceJunctionUnique(transData.TransName, geneReadData.getTranscripts(), prevExon.ExonEnd, exon.ExonStart);

                if(sjUnique)
                    ++uniqueSpliceJunctions;

                if(sjReads[TRANS_COUNT] > 0)
                {
                    ++spliceJunctionsSupported;

                    if(sjUnique)
                        ++uniqueSpliceJunctionsSupported;
                }
            }
        }

        double uniqueBaseAvgDepth = uniqueExonicBases > 0 ? uniqueExonicBaseTotalDepth / (double)uniqueExonicBases : 0;

        int[][] supportingFragments = geneCollection.getTranscriptReadCount(transData.TransId);

        double effectiveLength = calcEffectiveLength(exonicBases, expRateFragmentLengths);

        TranscriptResult results = ImmutableTranscriptResult.builder()
                .trans(transData)
                .exonsFound(exonsFound)
                .spliceJunctionsSupported(spliceJunctionsSupported)
                .uniqueSpliceJunctions(uniqueSpliceJunctions)
                .uniqueSpliceJunctionsSupported(uniqueSpliceJunctionsSupported)
                .spliceJunctionFragments(supportingFragments[FragmentMatchType.typeAsInt(FragmentMatchType.SPLICED)][TRANS_COUNT])
                .spliceJunctionUniqueFragments(supportingFragments[FragmentMatchType.typeAsInt(FragmentMatchType.SPLICED)][UNIQUE_TRANS_COUNT])
                .shortSupportingFragments(supportingFragments[FragmentMatchType.typeAsInt(FragmentMatchType.SHORT)][TRANS_COUNT])
                .shortUniqueFragments(supportingFragments[FragmentMatchType.typeAsInt(FragmentMatchType.SHORT)][UNIQUE_TRANS_COUNT])
                .longSupportingFragments(supportingFragments[FragmentMatchType.typeAsInt(FragmentMatchType.LONG)][TRANS_COUNT])
                .longUniqueFragments(supportingFragments[FragmentMatchType.typeAsInt(FragmentMatchType.LONG)][UNIQUE_TRANS_COUNT])
                .exonicBases(exonicBases)
                .exonicBaseCoverage(exonicBaseCoverage)
                .uniqueBases(uniqueExonicBases)
                .uniqueBaseCoverage(uniqueExonicBaseCoverage)
                .uniqueBaseAvgDepth(uniqueBaseAvgDepth)
                .effectiveLength(effectiveLength)
                .build();

        results.setFitAllocation(0);
        results.setTPM(0);

        return results;
    }

    public void setFitAllocation(double alloc) { mFitAllocation = alloc; }
    public void setTPM(double tpm) { mTPM = tpm; }
    public double getFitAllocation() { return mFitAllocation; }

    public static double calcEffectiveLength(int transLength, final List<int[]> fragmentLengthData)
    {
        if(fragmentLengthData.isEmpty())
            return transLength;

        long flFrequencyTotal = 0;
        long flBasesTotal = 0;

        for(final int[] flData : fragmentLengthData)
        {
            int fragLength = flData[FL_LENGTH];
            int fragFrequency = flData[FL_FREQUENCY];

            if(fragLength >= transLength)
                continue;

            long possibleBases = transLength - fragLength;
            flFrequencyTotal += fragFrequency;
            flBasesTotal += possibleBases * fragFrequency;
        }

        if(flFrequencyTotal == 0)
            return 0;

        return flBasesTotal / (double)flFrequencyTotal;
    }

    public double fragmentsPerKb()
    {
        return effectiveLength() > 0 ? mFitAllocation / (effectiveLength() / 1000.0) : 0;
    }

    private static boolean isSpliceJunctionUnique(final String transId, final List<TranscriptData> transDataList, long exonEnd, long exonStart)
    {
        for(TranscriptData transData : transDataList)
        {
            if (transData.TransName.equals(transId))
                continue;

            for (int i = 1; i < transData.exons().size(); ++i)
            {
                if(transData.exons().get(i-1).ExonEnd == exonEnd && transData.exons().get(i).ExonStart == exonStart)
                    return false;
            }
        }

        return true;
    }

    public static final String FLD_FIT_ALLOCATION = "FitAllocation";
    public static final String FLD_EFFECTIVE_LENGTH = "EffectiveLength";
    public static final String FLD_TPM = "TPM";

    public static String csvHeader()
    {
        return new StringJoiner(DELIMITER)
                .add(FLD_GENE_ID)
                .add(FLD_GENE_NAME)
                .add(FLD_TRANS_ID)
                .add(FLD_TRANS_NAME)
                .add("Canonical").add("ExonCount")
                .add(FLD_EFFECTIVE_LENGTH)
                .add(FLD_FIT_ALLOCATION)
                .add(FLD_TPM)
                .add("ExonsMatched").add("ExonicBases").add("ExonicCoverage")
                .add("UniqueBases").add("UniqueBaseCoverage").add("UniqueBaseAvgDepth")
                .add("SpliceJuncSupported").add("UniqueSpliceJunc").add("UniqueSpliceJuncSupported")
                .add("ShortFragments").add("ShortUniqueFragments").add("LongFragments").add("LongUniqueFragments")
                .add("SpliceJuncFragments").add("UniqueSpliceJuncFragments").toString();
    }

    public String toCsv(final EnsemblGeneData geneData)
    {
        return new StringJoiner(DELIMITER)
                .add(geneData.GeneId)
                .add(geneData.GeneName)
                .add(String.valueOf(trans().TransId))
                .add(trans().TransName)
                .add(String.valueOf(trans().IsCanonical))
                .add(String.valueOf(trans().exons().size()))
                .add(String.format("%.0f", effectiveLength()))
                .add(String.format("%.1f", mFitAllocation))
                .add(String.format("%6.3e", mTPM))
                .add(String.valueOf(exonsFound()))
                .add(String.valueOf(exonicBases()))
                .add(String.valueOf(exonicBaseCoverage()))
                .add(String.valueOf(uniqueBases()))
                .add(String.valueOf(uniqueBaseCoverage()))
                .add(String.valueOf(uniqueBaseAvgDepth()))
                .add(String.valueOf(spliceJunctionsSupported()))
                .add(String.valueOf(uniqueSpliceJunctions()))
                .add(String.valueOf(uniqueSpliceJunctionsSupported()))
                .add(String.valueOf(shortSupportingFragments()))
                .add(String.valueOf(shortUniqueFragments()))
                .add(String.valueOf(longSupportingFragments()))
                .add(String.valueOf(longUniqueFragments()))
                .add(String.valueOf(spliceJunctionFragments()))
                .add(String.valueOf(spliceJunctionUniqueFragments()))
                .toString();
    }
}
